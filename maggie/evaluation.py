import pandas as pd
import numpy as np
from .tasks import quality_control as qc 
from os.path import join, exists
from .utils import parse_read_counts, parse_quality_control_table


def construct_report(manifest, add_qc=True, add_cp=True, add_abundance=True, read_counts=None):
    rprt = list()
    for e in manifest.task_out_dir:
        if exists(e):
            rprt.append(pd.read_csv(e,'output/binning_table.csv'))
        else:
            print('Missing binning table: ', e)
    rprt = pd.concat(rprt)
    if add_qc:
        qctbls = list()
        for i in range(len(manifest)):
            t = manifest.iloc[i, :]
            if exists(join(t.task_out_dir, 'output/qc_table.csv')):
                qctbl = parse_quality_control_table(join(t.task_out_dir, 'output/qc_table.csv'))
                qc_measures = [m for q in t.qc_tools for m in getattr(qc, 'QCTool'+q).measures]
                qc_measures += ['bin', 'task_name']
                qctbls.append(qctbl[qc_measures])
            else:
                print('Missing quality table: ', e)
        qctbls = pd.concat(qctbls)
        float_cols = qctbls.select_dtypes(include=['float']).columns
        qctbls[float_cols] = qctbls[float_cols].astype(np.float32)
        rprt = pd.merge(rprt, qctbls, left_on=['task_name', 'representative_bin'], right_on=['task_name', 'bin'],how='left')
    if add_cp:
        assert Exception # This needs fixing
        asm = manifest.assembler.iloc[0] 
        cptbls = pd.read_parquet(
            join(utls.PNP_CP_CACHE_MH if asm == 'megahit' else utls.PNP_CP_CACHE_MS, 'rcvgnm_contig_properties.parquet')
        )
        # make float16 - don't need much accuracy, most values in 0-1
        float_cols = cptbls.select_dtypes(include=['float']).columns
        cptbls[float_cols] = cptbls[float_cols].astype(np.float32)
        # drop to relevant samples
        cptbls = cptbls[cptbls['sample'].isin(manifest.target_sample)]
        rprt = pd.merge(rprt,cptbls.drop(['sample','contig', 'contig_length', 'is_isolate', 'genome_length'],axis=1),on=['key','genome'])

    if add_abundance:
        assert exists(read_counts)
        read_counts = parse_read_counts(read_counts)
        target_samples = manifest.target_sample.unique()
        simulation_dir = manifest.iloc[0,'simulation_dir']
        abundance_map = pd.concat(
            [
                pd.read_csv(join(simulation_dir, f'{s}_metagenome_spec.csv')) for s in target_samples
            ]
        )[['Sample_file', 'genome', 'StrainAbund']].rename({'Sample_file':'sample'}, axis=1)
        abundance_map['genome_n_reads'] = abundance_map.apply(
            lambda x: read_counts[x['sample']] * x.StrainAbund, axis=1
        ).astype(int)
        abundance_map = abundance_map.set_index(['sample', 'genome'])
        rprt['genome_abundance'] = pd.Series(zip(rprt['sample'], rprt['genome'])).map(abundance_map.StrainAbund)
        rprt['genome_n_reads'] = pd.Series(zip(rprt['sample'], rprt['genome'])).map(abundance_map.genome_n_reads)
    rprt.index = range(len(rprt))
    return rprt

def construct_genome_metrics(rprt):
    float_cols = rprt.select_dtypes(include=['float16']).columns
    rprt[float_cols] = rprt[float_cols].astype(np.float64)
    metrics = compute_Fscore_metrics(rprt, ['task_name', 'genome'])
    return metrics


def add_report_data(rprt, gm, manifest, qctools, is_wrapper):

    # there are per-genome measurements or variables that are important for downstream evaluations
    all_qc_measures = list()
    for q in qctools:
        all_qc_meausres += getattr(qc, 'QCTool'+q).measures

    gm.drop(
        all_qc_measures + [
            'genome_length', 'genome_abundance', 'genome_n_reads', 'is_isolate',
            'assembler', 'binning_mode', 'sample', 'assembler_options', 'binner_options', 
            'refiner', 'refiner_options', 'is_refiner', 'summary_name'
        ], 
        errors='ignore', axis=1, inplace=True
    )
    # we need to look at only TP and FN records. Here the genome information, length, abundance, etc
    # and bin information, refer to the 'genome'. FP records need to be avoided, as the genome information / length etc
    # are about the genome from which the contig belongs to, which is not the genome in 'genome' - this contig is a false positive for that genome
    # i.e it shouldn't be there. Perhaps a better way to make the report is to make this consistent, relative to the genome studied.
    rprt = rprt[rprt.TP | rprt.FN][
        all_qc_measures + [
            'genome','task_name','genome_length', 'genome_abundance', 'genome_n_reads', 'is_isolate'
        ]].drop_duplicates()
    gm.reset_index(drop=True, inplace=True)
    gm = pd.merge(gm,rprt, left_on=['task_name', 'genome'], right_on=['task_name', 'genome'])

    manifest.set_index('task_name', inplace=True)
    gm['binning_mode'] = gm.task_name.map(manifest.binning_mode)
    gm['assembler'] = gm.task_name.map(manifest.assembler)
    gm['sample'] = gm.task_name.map(manifest.target_sample)
    gm['assembler_options'] = gm.task_name.map(manifest.assembler_options)
    gm['binner'] = manifest.apply(
        lambda x: x.binner if not x.is_refiner else f'{x.refiner}({x.summary_name})'
    )
    gm['binner_options'] =  gm.task_name.map(manifest.apply(lambda x: x.binner_options if not x.is_refiner else x.refiner_options))
    return gm
    
def compute_Fscore_metrics(bt, groupby, property=None, coverage_based=True, selected_metric=None):
    metrics = [
        bt.groupby(groupby).progress_apply(lambda x: f(x,w)).reset_index().rename({0:'value'},axis=1).assign(metric='cov_' + n if coverage_based else n)
        for f, w, n in (
            [
                e for e in [(precision, coverage_based, 'pr'), (recall, coverage_based, 'rc'), (fscore, coverage_based, 'fs')] 
                if (selected_metric is None or e[-1] == selected_metric)
            ]
        )
    ]
    metrics = pd.concat(metrics)
    if property:
        metrics['property'] = property
        metrics['n_bases'] = bt.groupby(groupby).apply(lambda x: x.contig_length.sum()).reset_index()[0]
        metrics['invpctl'] = metrics[property]
        metrics['pctl'] = metrics[property.replace('invpctl', 'pctl')]
        metrics.drop(property,axis=1,inplace=True)
    return metrics

def recall(x, cov_based):
    if cov_based:
        # x.TP == True if the contig is in the bin representing the genome.
        # i.e, the entire contig is a TP. But, we need to go through all the 
        # TP-contigs, can calculate the non-overlapping basepairs. This is the bp-tp.
        tp_contigs = x[x.TP]
        if tp_contigs.empty:
            return 0
        tp = x[x.TP].groupby('ref').apply(calculate_coverage).sum()

        # the false negative bp's are every part of the genome that's not covered
        # i.e fn genome len - tp. Since recall is tp/(tp+fn), we return tp/genome len
        return tp / x.iloc[0]['genome_length']
    else:
        # If its not coverage based, we're going to scores relative to the 
        # assembl-able, mappable portion of the metagenome. We call a FNs bp from 
        # contigs that are not assigned to the representative bin, rather than
        # the part of the genome not covered. We also don't do TP as coverage based.  
        denom = (x.TP * x.contig_length).sum() + (x.FN * x.contig_length).sum()
        if denom == 0:
            return 0
        return (x.TP * x.contig_length).sum() / denom

def precision(x, cov_based):
    if cov_based:
        tp_contigs = x[x.TP]
        if tp_contigs.empty:
            return 0
        tp = x[x.TP].groupby('ref').apply(calculate_coverage).sum()
        denom = tp + (x.FP * x.contig_length).sum()
        if denom == 0:
            return 1
        return tp / denom
    else:
        denom = (x.TP * x.contig_length).sum() + (x.FP * x.contig_length).sum()
        if denom == 0:
            return 1
        return (x.TP * x.contig_length).sum() / denom

def fscore(x, cov_based):
    p = precision(x, cov_based)
    r = recall(x, cov_based)
    denom = p+r
    if denom == 0:
        return 0
    return 2*p*r/(p+r)

def calculate_coverage(intervals):
    s = intervals.sort_values("ref_start")[["ref_start", "ref_end"]].to_numpy()
    total = 0
    cur_s, cur_e = s[0]
    cur_s, cur_e = (cur_s, cur_e) if cur_s < cur_e else (cur_e, cur_s)
    for start, end in s[1:]:
        start, end = (start, end) if start < end else (end, start)
        if start > cur_e:
            total += cur_e - cur_s
            cur_s, cur_e = start, end
        else:
            cur_e = max(cur_e, end)
    return total + (cur_e - cur_s)

def recoverable_genome_set(score_files, min_rc, min_pr):
    scores = parse_genome_measurement_files(score_files, min_cov=0)
    scores_rec = scores.groupby(['binning_mode', 'binner', 'assembler', 'sample', 'genome']).filter(
        lambda x: (x[x.metric == 'cov_rc'].value.iloc[0] >= min_rc) and (x[x.metric == 'cov_pr'].value.iloc[0] >= min_pr)
    )
    recovered_genomes = scores_rec[['sample', 'genome']].drop_duplicates()
    scores = scores.groupby(['binning_mode', 'binner', 'assembler'],group_keys=False).apply(
        lambda x: pd.merge(recovered_genomes, x, how='left', on=['sample', 'genome'])
    ).reset_index(drop=True)
    scores = scores[scores.metric.notna()]
    return scores, recovered_genomes


def parse_genome_measurement_files(score_files, min_cov, metric=None, isolate_only=True):
    score  = list()
    for sf in score_files:
        scr = pd.read_parquet(sf)
        if metric:
            scr = scr[scr.metric == metric] 
        score.append(scr)
    score = pd.concat(score)
    print('pre filtering', score.shape[0])
    # filter for genomes above min coverage
    score['genome_coverage'] = (score['genome_n_reads'] * 250) / score.genome_length
    score = score[score.genome_coverage >= min_cov]
    if isolate_only:
        score = score[score.is_isolate]
    print('post filtering', score.shape[0])
    score = score.reset_index(drop=True)
    return score