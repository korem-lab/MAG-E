import pandas as pd
import numpy as np
import glob
from .tasks import quality_control as qc, get_binner_coordinate, binner_coordinate
from os.path import join, exists, basename
from .utils import parse_read_counts, parse_quality_control_table, get_contig_name, remove_fasta_ext


def construct_binning_table(t):
    btfl = f'{t.binner_dir}/output/binning_table.pqt'
    if exists(btfl):
        return pd.read_parquet(btfl)
    bins = glob.glob(join(t.binner_dir, 'outputs/bins/*.fasta'))

    ground_truth = pd.read_csv(join(t.assembler_dir, f'ground_truth_table.csv'))
    records = list()
    for bin in bins:
        with open(bin) as f:
            contigs = [get_contig_name(e) for e in f if e.startswith('>')]
            records += [
                (*get_binner_coordinate(t), t.target, f'{t.target}-{e}', remove_fasta_ext(basename(bin)), e) 
                for e in contigs
            ]
    # filter for contigs that have been measured
    binned_contigs = pd.DataFrame(records, columns=[*list(binner_coordinate), 'target', 'key', 'bin', 'contig'])
    binned_contigs = binned_contigs.loc[binned_contigs.key.isin(ground_truth.key), :]
    binned_contigs.drop(['target', 'contig'],axis=1,inplace=True) # which contigs are in which bins

    # build binning table from ground truth
    binning_table = ground_truth[
        ['contig','contig_length','genome', 'sample', 'key','is_isolate', 
        'ref', 'ref_start', 'ref_end', 'genome_length']
    ].drop_duplicates()

    # Find the representative bin of each genome
    prod = pd.merge(binning_table, binned_contigs, on='key', how='outer')
    prod['representative_bin'] = prod.genome.map(prod.groupby('genome').apply(genome2bin_map))

    # add the representative bin to the binning table
    reprmap = prod[['genome', 'representative_bin']].drop_duplicates()
    reprmap = reprmap.set_index('genome', inplace=True)
    binning_table['representative_bin'] = binning_table.genome.map(reprmap.representative_bin)

    # compute contig-level TP and FN
    prod['contig_in_repr_bin'] = prod.bin == prod.representative_bin
    # In the product, its possible that g1, c1 ends up in three bins (if the binner allows it)
    # therefore, we groupby by g,c and ask if any of the bins the contig ends up in is the representative
    # of the genome. 
    tp = prod.groupby(['genome','key']).apply(lambda x: x.contig_in_repr_bin.any()).reset_index().rename({0:'ctgTP'},axis=1)
    fn = prod.groupby(['genome','key']).apply(lambda x: not x.contig_in_repr_bin.any()).reset_index().rename({0:'ctgFN'},axis=1)
    binning_table = binning_table.merge(tp,on=['genome','key'])
    binning_table = binning_table.merge(fn,on=['genome','key'])
    binning_table['ctgFP'] = False

    # ctgFP: go through each bin, and check the contig in the bin derives from the genome the bin represents.
    # If so, not a FP, if not its a FP. Since this is requires bin groupings, and a scan through ground truth
    # this requires a separate operation.
    fp = prod.groupby('bin').apply(compute_contig_fp).reset_index(drop=True)
    fp = fp[['genome','contig','contig_length','sample','key', 'ref', 'ref_start', 'ref_end', 'genome_length']].drop_duplicates()
    fp['ctgFP'] = True
    fp['ctgTP'] = False 
    fp['ctgFN'] = False
    # add isolate information to the genome
    isomap = ground_truth[['genome','is_isolate']].drop_duplicates()
    isomap.set_index('genome',inplace=True)
    fp['is_isolate'] = fp.genome.map(isomap.is_isolate)
    binning_table = pd.concat([binning_table, fp])

    # Write and return
    binning_table.to_parquet(btfl)
    return binning_table

def construct_quality_control_tables(t):
    tables = list()
    for tool in t.qctools:
        tool = getattr(qc, f'{tool}QCTool')
        tables.append(tool.to_table(t.binner_dir))
    return tables

def genome2bin_map(bt, use_genome_coverage=True):
    # for each contig in the bin, get the set of genomes it maps to
    if use_genome_coverage:
        # calculate the genome coverage along each chunk of the reference genome
        # that is found in the bin
        bin_sum = bt.groupby(['bin', 'ref'], dropna=False).apply(calculate_coverage).reset_index()
        # sum up over the bins
        bin_sum = bin_sum.groupby('bin', dropna=False)[0].sum()
    else:
        bin_sum = bt.groupby('bin',dropna=False).contig_length.sum()
    if bin_sum.index.notna().any():
        return bin_sum[bin_sum.index.notna()].idxmax()
    else:
        return pd.NA

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

def compute_contig_fp(df):
    # select for contigs from genomes that are represented by the bin.
    # A contig could be in a bin where the contig comes from a genome that the bin does not represent
    # and it should be counted as a false positive for the genome the bin does represent.
    filt=df.groupby('genome').filter(lambda x: all(x['bin'] == x['representative_bin']))
    if filt.empty:
        return pd.DataFrame(columns=['genome','contig','contig_length','sample','key'])
    fps=filt.groupby('genome').apply(lambda x: df[~df.key.isin(x.key)])
    fps.drop('genome',axis=1,inplace=True)
    fps=fps.reset_index()
    fps.drop('level_1',axis=1,inplace=True)
    return fps 

def construct_genome_metrics(binning_table, qctbls, read_count, ecodb, spec):

    # Compute the per-genome precision, recall, fscore metrics
    metrics = compute_Fscore_metrics(binning_table, groupby=['genome'])

    # add quality control values for each bin
    for tbl in qctbls:
        # why left: each genome assigned a representative bin, we care about that. 
        # We do not care (per se) about bins not assigned to be representative, hence left.
        metrics = pd.merge(metrics, tbl, left_on='representative_bin',right_on='bin', how='left')

    # add other genome properties
    isolates = set(ecodb.genome[ecodb.GenomeType == 'Isolate'])
    metrics['is_isolate'] = metrics.genome.isin(isolates)
    metrics['genome_length'] = metrics.genome.map(
        ecodb[['genome','Length']].set_index('genome').Length
    )
    metrics['genome_abundance'] = metrics.genome.map(
        spec[['genome', 'StrainAbund']].set_index('genome').StrainAbund
    )
    metrics['genome_n_reads'] = metrics.genome_abundance * read_count
    metrics['genome_coverage'] = (metrics.genome_n_reads * 250) / metrics.genome_length
    # E.G multi-strain, GC, TRNA etc etc
    return metrics


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

def calc_per_genome_metrics(self):
    # Construct the binning table for each binning tasks
    manifest = Manifest(self.config.project_base)
    tu.construct_binning_tables(manifest)

    # Build the reports. There is a report for each assembler, binner pair. 
    # A report contains the relevant information (genome, abundance, assigned bin, bin quality control)
    # for all (assembler, binner)-combo bin tasks. These reports are used to calculate the genome metrics.
    os.makedirs(self.config.evaluation_dir, exist_ok=True)
    bnab = manifest[['binner', 'assembler', 'is_refiner']].dropna().drop_duplicates().values
    rfab = manifest[['refiner', 'assembler', 'is_refiner']].dropna().drop_duplicates().values
    combos = np.vstack([bnab, rfab])
    for i in range(combos.shape[0]):
        bn, ab, is_refiner = combos[i]
        if is_refiner:
            _mnfst = manifest.loc[(manifest.refiner == bn) & (manifest.assembler == ab)]
        else:
            _mnfst = manifest.loc[(manifest.binner == bn) & (manifest.assembler == ab)]
        rprt = ev.construct_report(_mnfst, add_cp=False, read_counts=self.config.read_counts)
        rprt.to_parquet(join(self.config.evaluation_dir, f'{ab}_{bn}_REPORT.parquet'))

    # Compute the per-genome metrics
    for i in range(combos.shape[0]):
        bn, ab, _ = combos[i]
        rprt = pd.read_parquet(join(self.config.evaluation_dir, f'{ab}_{bn}_REPORT.parquet'))
        gm = ev.construct_genome_metrics(rprt)
        # these got dropped when constructing the genome metrics, all were doing here is adding them back
        gm = ev.add_report_data(rprt, gm, manifest.copy())
        gm.to_parquet(join(self.config.evaluation_dir, f'{ab}_{bn}.genome_metrics.parquet'))

def calc_contig_level_metrics(self, precision, recall, drop_raw_property=True):
    manifest = Manifest(self.config.project_base)
    # Add ground truth genome origin to the contig properties
    asm_tasks = manifest.loc[['assembler', 'target_sample', 'assembly_options', 'assembler_dir']].drop_duplicates()
    asm_tasks = asm_tasks[asm_tasks.assembly_options == 'default']
    cp_df_dataset = list()
    for i in range(len(asm_tasks)):
        t = asm_tasks.iloc[i,:]
        cp_df = parse_contig_properties(join(t.assembler_dir, f'{t.target_sample}_contig_properties.csv'))
        gt_df = pd.read_csv(join(t.assembler_dir, f'{t.target_sample}_gt_table.csv'))
        cp_df = pd.merge(cp_df, gt_df[['genome', 'key']], on='key')
        cp_df_dataset.append(cp_df)
    cp_df_dataset = pd.concat(cp_df_dataset)
    cp_df_dataset.reset_index(drop=True, inplace=True)

    # collect the genome metrics over the default tasks
    dflt_bin_tasks = manifest[
        (manifest.binner_options == 'default') & (manifest.assembly_options == 'default')
    ]
    gnm_metrics = ev.parse_genome_measurement_files(
        join(self.config.evaluation_dir, f'*.genome_metrics.parquet')
    )
    gnm_metrics = gnm_metrics.loc[gnm_metrics.task_hash.isin(dflt_bin_tasks.task_hash)]

    # get the recoverable set and filter to just those 
    _, rcvgnms = ev.recoverable_genome_set(gnm_metrics, recall, precision)
    cp_df_dataset = pd.merge(cp_df_dataset, rcvgnms, on=['sample', 'genome'])

    # iterate over continuous propertis and make a percentile form
    for prop in [e for e in cp_df_dataset.columns if e.startswith('prop_c')]:
        cp_df_dataset = pd.merge(
            cp_df_dataset, ev.to_percentiles(cp_df_dataset[[prop, 'contig_length']]),
            left_index=True, right_index=True, how='left'
        )
        if drop_raw_property:
            cp_df_dataset.drop(prop, axis=1, inplace=True)
    cp_df_dataset.to_parquet(self.config.evaluation_dir, 'processed_contig_properties.parquet')

    # break the eval up by binner and assembler combos
    bnab = manifest[['binner', 'assembler', 'is_refiner']].dropna().drop_duplicates().values
    rfab = manifest[['refiner', 'assembler', 'is_refiner']].dropna().drop_duplicates().values
    combos = np.vstack([bnab, rfab])
    for i in range(combos.shape[0]):
        bn, ab, is_refiner = combos[i]
        if is_refiner:
            _mnfst = manifest.loc[(manifest.refiner == bn) & (manifest.assembler == ab)]
        else:
            _mnfst = manifest.loc[(manifest.binner == bn) & (manifest.assembler == ab)]
        rprt = ev.construct_report(_mnfst, add_cp=True, add_abundance=False)
        rprt.to_parquet(join(self.config.evaluation_dir, f'{ab}_{bn}_CPREPORT.parquet'))
        ev.construct_contig_property_metrics(
            join(self.config.evaluation_dir, f'{ab}_{bn}_CPREPORT.parquet'), rcvgnms, manifest
        )

def evaluate_pipelines(self, precision, recall, plots=False):
    # Write the set of recoverable genome scores to a csv for analysis in R
    os.makedirs(join(self.config.evaluation_dir, 'LMM'), exist_ok=True)
    gnm_metrics = ev.parse_genome_measurement_files(
        glob.glob(join(self.config.evaluation_dir, '*.genome_metrics.parquet')), min_cov=0
    )
    scores, _ = ev.recoverable_genome_set(gnm_metrics, 0.2, 0.2)
    scores = scores[['binning_mode', 'binner', 'assembler', 'genome', 'is_refiner', 'sample', 'metric', 'value']]
    scores.to_csv(join(self.config.evaluation_dir, 'LMM/per_genome_metrics.csv'), index=None)

    # Run the LMM and write output to disk
    run_R_script('LMM_optimal_mag_pipeline', join(self.config.evaluation_dir, 'LMM'), 'true')
    ## plot data if needed
    #if plots:
    #    for by in ['binner', 'binning-mode', 'assembler']:
    #        for metric in ['fscore', 'precision', 'recall']:
    #            order = ['MaxBin2', 'VAMB', 'CONCOCT', 'METABAT2', 'SemiBin2', 'COMEBin']
    #            ax = plot_pipeline_performance_model(
    #                pd.read_csv(join(self.config.evaluation_dir, 'LMM', f'all_pipelines_{metric}.csv')), by=by, y_name=metric, order=order
    #            )
    #            plot_unlabelled_version(ax, join(utls.MANU_FIGS, 'Fig2', 'MAG_pipeline_LMM_fscore_binner'),ylower=-0.05, yupper=1.3,show=False)