import pandas as pd
import numpy as np
from .tasks import quality_control as qc 
from os.path import join, exists, basename
from .utils import parse_read_counts, parse_quality_control_table


def construct_report(manifest, add_qc=True, add_cp=True, add_abundance=True, read_counts=None, contig_props=None):
    rprt = list()
    for e in manifest.binner_dir:
        if exists(join(e, 'output/binning_table.csv')):
            rprt.append(pd.read_csv(join(e,'output/binning_table.csv'), index_col=0))
        else:
            print('Missing binning table: ', e)
    rprt = pd.concat(rprt)
    if add_qc:
        qctbls = list()
        for i in range(len(manifest)):
            t = manifest.iloc[i, :]
            if exists(join(t.binner_dir, 'output/qc_table.csv')):
                qctbl = parse_quality_control_table(join(t.binner_dir, 'output/qc_table.csv'))
                qc_measures = [m for q in t.qctools for m in getattr(qc, q+'QCTool')().qc_measures]
                qc_measures += ['bin', 'task_hash']
                qctbls.append(qctbl[qc_measures])
            else:
                print('Missing quality table: ', e)
        qctbls = pd.concat(qctbls)
        float_cols = qctbls.select_dtypes(include=['float']).columns
        qctbls[float_cols] = qctbls[float_cols].astype(np.float32)
        rprt = pd.merge(rprt, qctbls, left_on=['task_hash', 'representative_bin'], right_on=['task_hash', 'bin'],how='left')
        assert rprt[rprt.TP & rprt.bin.isna()].empty
    if add_cp:
        cptbls = pd.read_parquet(contig_props)
        # make float16 - don't need much accuracy, most values in 0-1
        float_cols = cptbls.select_dtypes(include=['float']).columns
        cptbls[float_cols] = cptbls[float_cols].astype(np.float32)
        rprt = pd.merge(rprt,cptbls.drop(['sample','contig', 'contig_length'],axis=1),on=['key','genome'])

    if add_abundance:
        assert exists(read_counts)
        read_counts = parse_read_counts(read_counts)
        target_samples = manifest.target_sample.unique()
        simulation_dir = manifest.iloc[0].simulation_dir
        abundance_map = pd.concat(
            [
                pd.read_csv(join(simulation_dir, f'{s}_metagenome_spec.csv')) for s in target_samples
            ]
        )[['Sample_file', 'genome', 'StrainAbund']].rename({'Sample_file':'sample'}, axis=1)
        abundance_map['genome_n_reads'] = abundance_map.apply(
            lambda x: read_counts.loc[basename(x['sample'])].item() * x.StrainAbund, axis=1
        ).astype(int)
        abundance_map['sample'] = abundance_map['sample'].apply(basename)
        abundance_map = abundance_map.set_index(['sample', 'genome'])
        rprt['genome_abundance'] = pd.Series(zip(rprt['sample'], rprt['genome'])).map(abundance_map.StrainAbund)
        rprt['genome_n_reads'] = pd.Series(zip(rprt['sample'], rprt['genome'])).map(abundance_map.genome_n_reads)
    rprt.index = range(len(rprt))
    return rprt

def construct_genome_metrics(rprt):
    float_cols = rprt.select_dtypes(include=['float16']).columns
    rprt[float_cols] = rprt[float_cols].astype(np.float64)
    metrics = compute_Fscore_metrics(rprt, ['task_hash', 'genome'])
    return metrics

def add_report_data(rprt, gm, manifest):

    # there are per-genome measurements or variables that are important for downstream evaluations
    all_qc_measures = list()
    qctools = list(set([x for l in manifest.qctools.drop_duplicates() for x in l]))
    for q in qctools:
        all_qc_measures += getattr(qc, q+'QCTool')().qc_measures

    gm.drop(
        all_qc_measures + [
            'genome_length', 'genome_abundance', 'genome_n_reads', 'is_isolate',
            'assembler', 'binning_mode', 'sample', 'assembler_options', 'binner_options',
            'refiner', 'refiner_options'
        ], 
        errors='ignore', axis=1, inplace=True
    )
    # we need to look at only TP and FN records. Here the genome information, length, abundance, etc
    # and bin information, refer to the 'genome'. FP records need to be avoided, as the genome information / length etc
    # are about the genome from which the contig belongs to, which is not the genome in 'genome' - this contig is a false positive for that genome
    # i.e it shouldn't be there. Perhaps a better way to make the report is to make this consistent, relative to the genome studied.
    rprt = rprt[rprt.TP | rprt.FN][
        all_qc_measures + [
            'genome','task_hash','genome_length', 'genome_abundance', 'genome_n_reads', 'is_isolate',
        ]].drop_duplicates()
    gm.reset_index(drop=True, inplace=True)
    gm = pd.merge(gm,rprt, left_on=['task_hash', 'genome'], right_on=['task_hash', 'genome'])

    manifest.set_index('task_hash', inplace=True)
    gm['binning_mode'] = gm.task_hash.map(manifest.binning_mode)
    gm['assembler'] = gm.task_hash.map(manifest.assembler)
    gm['sample'] = gm.task_hash.map(manifest.target_sample)
    gm['assembler_options'] = gm.task_hash.map(manifest.assembler_options)
    gm['is_refiner'] = gm.task_hash.map(manifest.is_refiner)
    gm['binner'] = gm.task_hash.map(manifest.apply(
        lambda x: f'{x.binner}{x.binner_summary_name}' if not x.is_refiner else f'{x.refiner}{x.refiner_summary_name}', axis=1
    ))
    gm['assembler'] = gm.task_hash.map(manifest.apply(lambda x: f'{x.assembler}{x.assembler_summary_name}', axis=1))
    gm['binner_options'] =  gm.task_hash.map(manifest.apply(lambda x: x.binner_options if not x.is_refiner else x.refiner_options, axis=1))
    return gm
    
def compute_Fscore_metrics(bt, groupby, property=None, coverage_based=True, selected_metric=None):
    metrics = [
        bt.groupby(groupby).apply(lambda x: f(x,w)).reset_index().rename({0:'value'},axis=1).assign(metric='cov_' + n if coverage_based else n)
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

def recoverable_genome_set(scores, min_rc, min_pr):
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
    # filter for genomes above min coverage
    score['genome_coverage'] = (score['genome_n_reads'] * 250) / score.genome_length
    score = score[score.genome_coverage >= min_cov]
    if isolate_only:
        score = score[score.is_isolate]
    score = score.reset_index(drop=True)
    return score

def to_percentiles(data):
    data = data.dropna()
    name = data.columns[0] # prop name
    data.sort_values(by=name, inplace=True)
    data['clen_csum'] = data.contig_length.cumsum()
    percentiles = (np.arange(start=1, stop=101) * data.contig_length.sum()/100).astype(int)
    data[f'{name}_pctl'] = data.clen_csum.apply(lambda x: np.searchsorted(percentiles, x)).astype(int)

    # downgrade the percentile if the majority (> 50% length) of a contig is present
    # in the lower percentile
    data[f'{name}_pctl'] = data.progress_apply(
        lambda x: x[f'{name}_pctl']
        if x[f'{name}_pctl'] == 0 else x[f'{name}_pctl'] - 1
        if ((x.clen_csum - percentiles[int(x[f'{name}_pctl'])-1]) < x.contig_length/2.0)
        else x[f'{name}_pctl'], axis=1
    ).astype(int)

    # compute the last element of each percentile .
    # -1 because we want the last element of each percentile, not the first of each
    invpctl = data.iloc[
        np.searchsorted(data[f'{name}_pctl'].values, np.arange(100),side='right')-1
    ][name].reset_index(drop=True)
    # map each percentile (0-99) to its inverse 
    data[f'{name}_invpctl'] = data[f'{name}_pctl'].map(invpctl)
    return data[[f'{name}_pctl', f'{name}_invpctl']]


def construct_contig_property_metrics(prq_file, rcvgnms, manifest):
    # read report 
    rprt = pd.read_parquet(prq_file)

    # calculate continuous property metrics
    rprt = pd.merge(rprt, rcvgnms, on=['sample', 'genome'])
    all_metrics = list()
    for p in [e for e in rprt.columns if 'invpctl' in e]:
        assert not any(np.isinf(rprt[p].values))
        metrics = compute_Fscore_metrics(rprt, ['task_hash', p, p.replace('invpctl', 'pctl')], property=p, selected_metric='rc', coverage_based=False)
        all_metrics.append(metrics)
    all_metrics = pd.concat(all_metrics)

    # calculate discrete property metrics
    rprt = pd.merge(
        manifest, bntsks[['binner','binning_mode','task_hash']].drop_duplicates(), on='task_hash'
    )

    all_discrete_metrics = list()
    for p in [e for e in rprt.columns if 'prop_d' in e]:
        levels = rprt[p].unique()
        scores = list()
        for level in levels:
            score_level = compute_Fscore_metrics(
                rprt[rprt.p == level], 
                groupby=['task_hash', 'binning_mode', 'binner', 'assembler', 'sample', 'genome'], 
                selected_metric='rc', coverage_based=False
            )
            score_level['level'] = level
            scores.append(score_level)
        scores = pd.concat(scores)
        scores['property'] = p
        all_discrete_metrics.append(scores)
    all_metrics = pd.concat([all_metrics, all_discrete_metrics])
    all_metrics = all_metrics[
        [
            'task_hash', 'value','metric','property','n_bases','invpctl','pctl', 'genome',
            'binning_mode', 'binner', 'assembler', 'sample', 'level'
        ]
    ]
    all_metrics.to_parquet(prq_file.replace('_CPREPORT.parquet', f'.cp_metrics.parquet'))