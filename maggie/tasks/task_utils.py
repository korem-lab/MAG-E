from os.path import join, exists, basename
import pandas as pd
import hashlib
from functools import reduce
from itertools import product
from collections import defaultdict
import json
from ..utils import parse_binning_mode_datasets, rm_dir, get_contig_name, remove_fasta_ext, flatten, run
from . import assembly as ab, binning as bn, quality_control as qc, mapping as mp


def get_summary_names(data, is_refiner=False):
    print(data)
    if is_refiner:
        summary_names = pd.DataFrame(
            {'tool':refiner, 'toolopt':ro, 'pipelines':pipelines} for refiner, ro, pipelines in data
        )
    else:
        summary_names = pd.DataFrame({'tool':t, 'toolopt':to} for (t,to) in data)
    summary_names['sname'] = summary_names.groupby('tool').cumcount().add(1).astype(str)
    summary_names['sname'] = summary_names['sname'].apply(lambda x: f'({x})')
    mask = summary_names.groupby('tool')['tool'].transform('count') == 1
    summary_names.loc[mask, 'sname'] = ''
    return summary_names

def id_gen():
    counter = 1
    def new_id():
        nonlocal counter
        val = counter
        counter += 1
        return val
    return defaultdict(new_id)

def make_ids(l):
    ids = dict()
    for k,v in l:
        if k in ids:
            ids[k][v] = len(ids[k])+1
        else:
            ids[k] = {v:1}
    return ids

def make_manifest(
    assemblers, mappers, binners, modes, refiners, qctools, simulation_dir, 
    assembly_cache, mapping_cache, bintask_dir, project_base
):
    # The same tool can be provided, but with different options. 
    # Mapping each (tool, option) pair to a numerical id gives a simple naming scheme
    aids = make_ids(assemblers)
    mids = make_ids(mappers)
    bids = make_ids(binners)
    id_to_string = lambda ids, t, o: f'({ids[t][o]})' if len(ids[t])>1 else ""

    # Enumerate the different directories 
    bintask_id = id_gen()
    asmcache_id = id_gen()
    mapcache_id = id_gen()
    tasks = list()
    for (ab, abo), (mp, mpo), (bn, bno), mode in product(assemblers, mappers, binners, modes):

        asm_num = asmcache_id[(ab, abo)]
        map_num = mapcache_id[(ab, abo, mp, mpo)]
        bin_task_num = bintask_id[(ab, abo, mp, mpo, bn, bno, mode)]
        spec = {
            'simulation_dir': simulation_dir,
            'assembler': ab, 'assembler_options': abo, 'binner': bn, 'binner_options': bno,
            'mapper': mp, 'mapper_options': mpo,
            'binning_mode': mode, 'qctools': qctools,
            'is_refiner': False
        }

        # get the summary names for the binning runs
        spec['assembler_summary_name'] = id_to_string(aids, ab, abo)
        spec['binner_summary_name'] = id_to_string(bids, bn, bno)
        spec['mapping_summary_name'] = id_to_string(mids, mp, mpo)

        datasets = parse_binning_mode_datasets(join(project_base, f'{mode}_datasets.csv'))
        for trgt, df in datasets.groupby('target'):
            spec['target'] = trgt
            spec['samples'] = [trgt] + sorted(list(set(df.dataset.to_list()) - {trgt}))
            spec['task_out_dir'] = join(bintask_dir, f'bin_task_{bin_task_num}', trgt)
            spec['asm_dir'] = join(assembly_cache, f'assembly_task_{asm_num}', trgt)
            spec['map_dir'] = join(mapping_cache, f'mapping_task_{map_num}', trgt)
            tasks.append(spec.copy())

    # Refiners will be considered just another binner.
    # i.e, it will the refiner and options will be listed under binner and binner options
    # however, the is_refiner will be true, and their will be a pipeline column
    for (rf, ro, ab, abo, mp, mpo, binner_set) in refiners:
        bin_task_num = bintask_id[(rf, ro, tuple(binner_set))]
        asm_num = asmcache_id((ab, abo))
        map_num = mapcache_id((ab, abo, mp, mpo))
        spec = {
            'simulation_dir': simulation_dir,
            'assembler': ab, 'assembler_options': abo, 'binner': rf, 'binner_options': ro,
            'mapper': mp, 'mapper_options': mpo, 'binner_set': binner_set, 'is_refiner':True,
            'qctools': qctools
        }
        spec['binner_summary_name'] = id_to_string(bids, rf, ro+' '.join(flatten(binner_set)))
        for trgt, df in datasets.groupby('target'):
            spec['target'] = trgt
            spec['samples'] = pd.NA
            spec['task_out_dir'] = join(bintask_dir, f'bin_task_{bin_task_num}', trgt)
            spec['asm_dir'] = join(assembly_cache, f'assembly_task_{asm_num}', trgt)
            tasks.append(spec.copy())
            spec['map_dir'] = join(mapping_cache, f'mapping_task_{map_num}', trgt)
    for spec in tasks:
        spec['task_name'] = hashlib.sha256(json.dumps(spec).encode()).hexdigest()
    tasks = pd.DataFrame(tasks)
    return tasks

def get_tasks(
        manifest, assembler, aopt, mapper, mopt, binner, bopt, binning_mode, binner_set, target
    ):
    flt = lambda x,y: manifest[x] == y if y is not None else pd.Series(True, index=manifest.index)
    if binner_set:
        binner_set = tuple(tuple(e.split(',')) for e in binner_set.split(':'))
    return manifest.loc[
        flt('target', target) & flt('assembler', assembler) & flt('assembler_options', aopt) &
        flt('mapper', mapper) & flt('mapper_options', mopt) & flt('binner', binner) & 
        flt('binner_options', bopt) & flt('binning_mode', binning_mode) &
        flt('binner_set', binner_set)
    ]

def get_task_directory(self, **kwargs):
    tsks = self.get_tasks(**kwargs)
    assembler, mapper, binner = kwargs['assembler'], kwargs['mapper'], kwargs['binner']
    binning_mode, target = kwargs['binning_mode'], kwargs['target']
    if assembler and mapper and binner and binning_mode and target:
        tsks[['task_out_dir']].drop_duplicates()
    elif assembler and mapper and target:
        tsks[['map_dir']].drop_duplicates()
    elif assembler and target:
        tsks[['asm_dir']].drop_duplicates()
    else:
        Exception("Invalid query.")
    if tsks.empty:
        raise Exception("No tasks fit the query.")
    elif len(tsks) > 1:
        raise Exception("More than task fits the query") 
    return tsks.iloc[0,0]

def run_assembly(task, threads=8, force=False, check=False):
    r1 = join(task.simulation_dir, f'{task.target}_R1.fastq.gz')
    r2 = join(task.simulation_dir, f'{task.target}_R2.fastq.gz')
    options = '' if task.assembler_options == 'default' else task.assembler_options
    assembler = getattr(ab, f'{task.assembler}Assembler')()
    kwargs = task.to_dict() | {'threads':threads, 'r1':r1, 'r2':r2,'options':options}
    if check:
        return assembler.main_done(**kwargs)
    if force or not assembler.main_done(**kwargs):
        rm_dir(task.asm_dir, remake=True)
        assembler.run_main(**kwargs)

def run_mapping(task, stage, map_sample, threads=8, force=False, check=False):
    assert map_sample in task.samples
    r1 = join(task.simulation_dir, f'{map_sample}_R1.fastq.gz')
    r2 = join(task.simulation_dir, f'{map_sample}_R2.fastq.gz')
    options = '' if task.mapper_options == 'default' else task.mapper_options
    mapper = getattr(mp, f'{task.mapper}Mapper')()
    kwargs = {'r1':r1, 'r2':r2, 'map_sample':map_sample, 'threads':threads, 'options':options}
    kwargs = task.to_dict() | kwargs
    if stage == 'prep' or stage == 'all':
        if check:
            return mapper.prep_done(**kwargs)
        if force or not mapper.prep_done(**kwargs):
            rm_dir(task.map_dir, remake=True)
            mapper.run_prep(**kwargs)
    if stage == 'main' or stage == 'all':
        if check:
            return mapper.run_main(**kwargs)
        if force or not mapper.main_done(**kwargs):
            mapper.run_main(**kwargs)

def run_binning(task, stage, threads=8, force=False, check=False):
    binner = getattr(bn, f'{task.binner}Binner')()
    kwargs = task.to_dict() | {'threads':threads}
    options = '' if task.binner_options == 'default' else task.binner_options
    if stage == 'prep' or stage == 'all':
        if check:
            return binner.prep_done(**kwargs)
        if force or not binner.prep_done(**kwargs):
            rm_dir(join(task.task_out_dir, 'input'),remake=True)
            binner.run_prep(**kwargs)
    if stage == 'main' or stage == 'all':
        if check:
            return binner.main_done(**kwargs)
        if force or not binner.main_done(**kwargs):
            rm_dir(join(task.task_out_dir, 'output'),remake=True)
            binner.run_main(**kwargs)

def run_quality_control(task, qctool, threads=8, force=False, check=False):
    qctool = getattr(qc, f'{qctool}QCTool')()
    kwargs = task.to_dict() | {'threads':threads}
    if check:
        qctool.main_done(**kwargs)
    if force or not qctool.main_done(**kwargs):
        qctool.run_main(**kwargs)

#def run_quality_control(manifest, force=False, threads=8):
#    for i in range(len(manifest)):
#        t = manifest.loc[i,:]
#        tables = list()
#        for q in t.qctools:
#            q = getattr(qc, q+'QCTool')()
#            if not q.has_bins(t.task_out_dir):
#                continue
#            if not q.done(t.task_out_dir):
#                q.run(**(t.to_dict() | {'threads':threads}))
#            tables.append(q.to_qctable(**(t.to_dict() | {'threads':threads})))
#        if tables:
#            qc_table = reduce(lambda left,right: pd.merge(left,right,on='bin'),tables)
#            qc_table['task_name'] = t.task_name
#            qc_table.to_csv(join(t.task_out_dir, 'output/qc_table.csv'), index=None)

#def run_mappings(manifest, stage='all', threads=8):
#    map_tasks = manifest[['target', 'simulation_dir', 'asm_dir', 'map_dir','mapper', 'mapper_options', 'samples']].drop_duplicates()
#    for i in range(len(map_tasks)):
#        t = map_tasks.iloc[i,:]
#        for s in t.samples:
#            run_mapping(t, stage, map_sample=s, threads=threads)

#def run_assemblies(manifest, threads=8):
#    asm_tasks = manifest[~manifest.is_refiner][['target', 'simulation_dir', 'asm_dir', 'assembler', 'assembler_options']].drop_duplicates()
#    print('Total assembly tasks: ', len(asm_tasks))
#    for i in range(len(asm_tasks)):
#        run_assembly(asm_tasks.iloc[i,:], threads)

#def run_binnings(manifest, run_only=None, force_prep=False, force_bin=False, threads=8):
#    # run prep
#    bin_tasks = manifest[~manifest.is_refiner]
#    if run_only is None or run_only == 'prep':
#        for i in range(len(bin_tasks)):
#            t = bin_tasks.iloc[i,:]
#            if force_prep or not task_done(t, 'prep', clear=True):
#                run_prep(t, threads)
#    # run bin
#    if run_only is None or run_only == 'bin':
#        for i in range(len(bin_tasks)):
#            t = bin_tasks.iloc[i,:]
#            if force_bin or not task_done(t, 'bin', clear=True):
#                run_bin(t, threads)
#
#def run_refine_tasks(manifest, run_only=None, force_refine=False, threads=8):
#    refine_tasks = manifest[manifest.is_refiner] 
#    bin_tasks = manifest[~manifest.is_refiner]
#    if run_only is None or run_only == 'refine':
#        for i in range(len(refine_tasks)):
#            t = refine_tasks.iloc[i, :]
#            if force_refine or not task_done(t, 'bin', clear=True):
#                ppln_bin_out = list()
#                for (ab, abo, bn, bno, mode) in t.pipelines:
#                    # select bin task
#                    bt = bin_tasks.loc[
#                        (bin_tasks.assembler == ab) & (bin_tasks.assembler_options == abo) & 
#                        (bin_tasks.binner == bn) & (bin_tasks.binner_options == bno) & 
#                        (bin_tasks.binning_mode == mode) & (bin_tasks.target == t.target)
#                    ]
#                    assert len(bt) == 1
#                    # get the path
#                    ppln_bin_out.append(bt.task_out_dir.item())
#                run_bin(t, threads, refine=True, ppln_bin_out=ppln_bin_out)


def construct_binning_tables(manifest):
    for i in range(len(manifest)):
        t = manifest.iloc[i,:]
        if t.is_refiner:
            b = getattr(bn, t.refiner+'Refiner')()
        else:
            b = getattr(bn, t.binner+'Binner')()
        bins = b.bins_as_fasta(join(t.task_out_dir, 'output/bins'))
        bint = compute_binning_table(t.task_name, t.target, bins, t.asm_dir)
        bint.to_csv(join(t.task_out_dir, 'output', 'binning_table.csv'))

def compute_binning_table(task_name, sample, bins, gt_dir):
    ground_truth = pd.read_csv(join(gt_dir, f'{sample}_gt_table.csv'))
    ground_truth = ground_truth.loc[ground_truth['sample'] == sample]
    records = list()
    for bin in bins:
        with open(bin) as f:
            contigs = [get_contig_name(e) for e in f if e.startswith('>')]
            records += [
                (task_name, sample, f'{sample}-{e}', remove_fasta_ext(basename(bin)), e) 
                for e in contigs
            ]
    # filter for contigs that have been measured
    binned_contigs = pd.DataFrame(records, columns=['task_name', 'sample', 'key', 'bin', 'contig'])
    binned_contigs = binned_contigs.loc[binned_contigs.key.isin(ground_truth.key), :]
    binned_contigs.drop(['sample', 'task_name', 'contig'],axis=1,inplace=True) # which contigs are in which bins

    # build binning table from ground truth
    binning_table = ground_truth[
        ['contig','contig_length','genome', 'sample', 'key','is_isolate', 'ref', 'ref_start', 'ref_end', 'genome_length']
    ].drop_duplicates()

    prod = pd.merge(binning_table, binned_contigs, on='key', how='outer')
    prod['representative_bin'] = prod.genome.map(prod.groupby('genome').apply(genome2bin_map))
    # compute TP
    prod['contig_in_repr_bin'] = prod.bin == prod.representative_bin
    tp = prod.groupby(['genome','key']).apply(lambda x: x.contig_in_repr_bin.any()).reset_index().rename({0:'TP'},axis=1)
    fn = prod.groupby(['genome','key']).apply(lambda x: not x.contig_in_repr_bin.any()).reset_index().rename({0:'FN'},axis=1)
    binning_table = binning_table.merge(tp,on=['genome','key'])
    binning_table = binning_table.merge(fn,on=['genome','key'])
    binning_table['FP'] = False
    # Computing FP is more involved.
    fp = prod.groupby('bin').apply(compute_fp).reset_index(drop=True)
    fp = fp[['genome','contig','contig_length','sample','key', 'ref', 'ref_start', 'ref_end', 'genome_length']].drop_duplicates()
    fp['FP'] = True
    fp['TP'] = False 
    fp['FN'] = False
    # add isolate information to the genome
    isomap = ground_truth[['genome','is_isolate']].drop_duplicates()
    isomap.index=isomap.genome
    isomap.drop('genome',axis=1,inplace=True)
    fp['is_isolate'] = fp.genome.map(isomap.is_isolate)
    binning_table = pd.concat([binning_table, fp])

    # add back the representative bin of each genom
    reprmap = prod[['genome', 'representative_bin']].drop_duplicates()
    reprmap.index = reprmap.genome
    reprmap.drop('genome',axis=1,inplace=True)
    binning_table['representative_bin'] = binning_table.genome.map(reprmap.representative_bin)
    binning_table['task_name'] = task_name
    # First compute FP for genomes that have a bin representation
    return binning_table

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

def compute_fp(df):
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


def compute_idxstats(bam, idxstats):
    run(f'samtools idxstats {bam} > {idxstats}')