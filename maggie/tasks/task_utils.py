from os.path import join, exists, basename
import pandas as pd
import hashlib
from subprocess import run, DEVNULL
from functools import reduce
from itertools import product
from collections import defaultdict
import json
from ..utils import parse_binning_mode_datasets, write_done_flag, rm_dir, get_contig_name, remove_fasta_ext
from . import assembly as ab, binning as bn, quality_control as qc


def get_summary_names(data, dtype):
    assert dtype in ['refiner', 'bn', 'ab']
    if dtype == 'refiner':
        summary_names = pd.DataFrame(
            {'refiner':refiner, 'ro':ro, 'pipelines':pipelines} for refiner, ro, pipelines in data
        )
    elif dtype == 'bn':
        summary_names = pd.DataFrame({'bn': bn, 'bno':bno} for (bn, bno) in data)
    else:
        summary_names = pd.DataFrame({'ab': ab, 'abo':abo} for (ab, abo) in data)

    
    summary_names['sname'] = summary_names.groupby(dtype).cumcount().add(1).astype(str)
    summary_names['sname'] = summary_names['sname'].apply(lambda x: f'({x})')
    mask = summary_names.groupby(dtype)[dtype].transform('count') == 1
    summary_names.loc[mask, 'sname'] = ''
    return summary_names

def make_key_map():
    counter = 1
    def new_id():
        nonlocal counter
        val = counter
        counter += 1
        return val
    return defaultdict(new_id)

def make_manifest(
    assemblers, binners, modes, refiners, qctools, simulation_dir, 
    assembly_cache, bintask_dir, project_base
):
    # Assert their unique
    def chck(x):
        assert len(x) == len(set(x)), Exception("List myst be unique")
    chck(assemblers)
    chck(binners)
    chck(modes)
    chck(refiners)

    # There may be multiple assembler, binner, and refiners run with different options
    # to give these a simple name, which can later be looked up, we construct tables
    # that map each option to an integer "summary name"
    rsmry = get_summary_names(refiners, 'refiner')
    bsmry = get_summary_names(binners, 'bn')
    asmry = get_summary_names(assemblers, 'ab')

    bintask_counter = make_key_map()
    assembly_counter = make_key_map()
    tasks = list()
    for (ab, abo), (bn, bno), mode in product(assemblers, binners, modes):

        # abo and bno are respectively the assembly options and binning options. 
        # For each option copbination, we make a separate directory
        asm_num = assembly_counter[(ab, abo)]
        bin_task_num = bintask_counter[(ab, abo, bn, bno, mode)]
        spec = {
            'simulation_dir': simulation_dir,
            'assembler': ab, 'assembler_options': abo, 'binner': bn, 'binner_options': bno,
            'binning_mode': mode, 'qctools': qctools,
            'is_refiner': False
        }

        # get the summary names for the binning runs
        spec['assembler_summary_name'] = asmry[(asmry.ab == ab) & (asmry.abo == abo)].sname.item()
        spec['binner_summary_name'] = bsmry[(bsmry.bn == bn) & (bsmry.bno == bno)].sname.item()

        datasets = parse_binning_mode_datasets(join(project_base, f'{mode}_datasets.csv'))
        for trgt, df in datasets.groupby('target_sample'):
            spec['target_sample'] = trgt
            spec['samples'] = [trgt] + sorted(list(set(df.dataset.to_list()) - {trgt}))
            spec['task_out_dir'] = join(bintask_dir, f'bin_task_{bin_task_num}', trgt)
            spec['asm_dir'] = join(assembly_cache, f'assembly_task_{asm_num}', trgt)
            tasks.append(spec.copy())


    for (refiner, ro, pipelines) in refiners:
        bin_task_num = bintask_counter[(refiner, ro, tuple(pipelines))]
        asm_num = assembly_counter[(pipelines[0][0], pipelines[0][1])]
        spec = {
            'qctools': qctools, 'simulation_dir': simulation_dir, 'is_refiner': True,
            'refiner': refiner, 'refiner_options': ro, 'pipelines': pipelines,
            'assembler': pipelines[0][0], 'assembler_options':pipelines[0][1]
        }
        spec['refiner_summary_name'] = rsmry.loc[(rsmry.refiner == refiner) & (rsmry.ro == ro) & (rsmry.pipelines == pipelines)].sname.item()
        datasets = parse_binning_mode_datasets(join(project_base, f'{mode}_datasets.csv'))
        for trgt, df in datasets.groupby('target_sample'):
            spec['target_sample'] = trgt
            spec['samples'] = [trgt] + sorted(list(set(df.dataset.to_list()) - {trgt}))
            spec['task_out_dir'] = join(bintask_dir, f'bin_task_{bin_task_num}', trgt)
            spec['asm_dir'] = join(assembly_cache, f'assembly_task_{asm_num}', trgt)
            tasks.append(spec.copy())
            
    for spec in tasks:
        spec['task_name'] = hashlib.sha256(json.dumps(spec).encode()).hexdigest()
    tasks = pd.DataFrame(tasks)
    return tasks

def run_assemblies(manifest, threads=8):

    # multiple bintasks can use the sample assembly so first,
    # gather all the unique assembly tasks from the manifest
    asm_tasks = manifest[~manifest.is_refiner][['target_sample', 'simulation_dir', 'asm_dir', 'assembler', 'assembler_options']].drop_duplicates()
    print('Total assembly tasks: ', len(asm_tasks))
    # run tasks
    for i in range(len(asm_tasks)):
        t = asm_tasks.iloc[i,:]
        r1 = join(t.simulation_dir, f'{t.target_sample}_R1.fastq.gz')
        r2 = join(t.simulation_dir, f'{t.target_sample}_R2.fastq.gz')
        options = '' if t.assembler_options == 'default' else t.assembler_options
        assembler = getattr(ab, f'{t.assembler}Assembler')()
        assembler.run_assembly(r1, r2, t.asm_dir, threads, options)
        assembler.clean_up(t.target_sample, t.asm_dir)
        assert assembler.assembly_done(t.target_sample, t.asm_dir)

def run_mapping(manifest, threads=8):
    # multiple bintasks can use the sample assembly so first,
    # gather all the unique assembly tasks from the manifest
    manifest = manifest[~manifest.is_refiner]
    asm_tasks = manifest[['target_sample', 'simulation_dir', 'asm_dir', 'assembler', 'assembler_options']].drop_duplicates()
    for i in range(len(asm_tasks)):
        t = asm_tasks.iloc[i,:]
        samples = set([s for ss in manifest[manifest.target_sample == t.target_sample]['samples'] for s in ss])
        for s in samples:
            r1 = join(t.simulation_dir, f'{s}_R1.fastq.gz')
            r2 = join(t.simulation_dir, f'{s}_R2.fastq.gz')
            contigs = join(t.asm_dir, f'{t.target_sample}.fasta')
            bowtie2_build(contigs, contigs.replace('.fasta', ''))
            bam = contigs.replace('.fasta', f'_{s}.bam')
            bowtie2(contigs.replace('.fasta',''),bam, r1, r2, threads)
            sort_bam(bam,threads=threads)
            index_bam(bam)

def run_bin_tasks(manifest, run_only=None, force_prep=False, force_bin=False, threads=8):
    # run prep
    bin_tasks = manifest[~manifest.is_refiner]
    if run_only is None or run_only == 'prep':
        for i in range(len(bin_tasks)):
            t = bin_tasks.iloc[i,:]
            if force_prep or not task_done(t, 'prep', clear=True):
                run_prep(t, threads)
    # run bin
    if run_only is None or run_only == 'bin':
        for i in range(len(bin_tasks)):
            t = bin_tasks.iloc[i,:]
            if force_bin or not task_done(t, 'bin', clear=True):
                run_bin(t, threads)

def run_refine_tasks(manifest, run_only=None, force_refine=False, threads=8):
    refine_tasks = manifest[manifest.is_refiner] 
    bin_tasks = manifest[~manifest.is_refiner]
    if run_only is None or run_only == 'refine':
        for i in range(len(refine_tasks)):
            t = refine_tasks.iloc[i, :]
            if force_refine or not task_done(t, 'bin', clear=True):
                ppln_bin_out = list()
                for (ab, abo, bn, bno, mode) in t.pipelines:
                    # select bin task
                    bt = bin_tasks.loc[
                        (bin_tasks.assembler == ab) & (bin_tasks.assembler_options == abo) & 
                        (bin_tasks.binner == bn) & (bin_tasks.binner_options == bno) & 
                        (bin_tasks.binning_mode == mode) & (bin_tasks.target_sample == t.target_sample)
                    ]
                    assert len(bt) == 1
                    # get the path
                    ppln_bin_out.append(bt.task_out_dir.item())
                run_bin(t, threads, refine=True, ppln_bin_out=ppln_bin_out)

def run_prep(spec, threads):
    print('run_prep:', spec.task_name, flush=True)
    b = getattr(bn, spec.binner + 'Binner')()
    b.run_prep(**(spec.to_dict() | {'threads':threads}))
    write_done_flag(spec.task_out_dir, name='prep')

def run_bin(spec, threads, refine=False, ppln_bin_out=None):
    print('run_bin:', spec.task_name,flush=True)
    b = getattr(bn, (spec.binner + 'Binner') if not refine else (spec.refiner + 'Refiner'))()
    spec.binner_options = '' if spec.binner_options == 'default' else spec.binner_options
    spec.refiner_options = '' if spec.refiner_options == 'default' else spec.refiner_options
    if refine:
        assert ppln_bin_out 
        aux = {'ppln_bin_out':ppln_bin_out, 'threads':threads}
    else:
        aux = {'threads':threads}
    b.run_binning(**(spec.to_dict() | aux))
    write_done_flag(spec.task_out_dir, name='bin')

def clear_prep(task_out_dir):
    rm_dir(join(task_out_dir, 'input'))
    rm_dir(join(task_out_dir, 'prep_DONE'))

def clear_bin(task_out_dir):
    rm_dir(join(task_out_dir, 'output/bins'))
    rm_dir(join(task_out_dir, 'bin_DONE'))

def run_quality_control(manifest, force=False, threads=8):
    for i in range(len(manifest)):
        t = manifest.loc[i,:]
        tables = list()
        for q in t.qctools:
            q = getattr(qc, q+'QCTool')()
            if not q.has_bins(t.task_out_dir):
                continue
            if not q.done(t.task_out_dir):
                q.run(**(t.to_dict() | {'threads':threads}))
            tables.append(q.to_qctable(**(t.to_dict() | {'threads':threads})))
        if tables:
            qc_table = reduce(lambda left,right: pd.merge(left,right,on='bin'),tables)
            qc_table['task_name'] = t.task_name
            qc_table.to_csv(join(t.task_out_dir, 'output/qc_table.csv'), index=None)

def construct_binning_tables(manifest):
    for i in range(len(manifest)):
        t = manifest.iloc[i,:]
        if t.is_refiner:
            b = getattr(bn, t.refiner+'Refiner')()
        else:
            b = getattr(bn, t.binner+'Binner')()
        bins = b.bins_as_fasta(join(t.task_out_dir, 'output/bins'))
        bint = compute_binning_table(t.task_name, t.target_sample, bins, t.asm_dir)
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

def task_done(bntsk, stage='all', report=False, only_not_done=False, clear=False):
    # TODO ADD QC STAGE, and BIN TABLE STAGE
    if stage == 'all' or stage == 'prep':
        task_out_dir = bntsk.task_out_dir
        b = getattr(bn, bntsk.binner+'Binner')()
        prep_done =  b.prep_done(**bntsk.to_dict())
        prep_done &= exists(join(task_out_dir, 'prep_DONE'))
        if not prep_done and report:
            print('bin prep not DONE', bntsk.task_name, flush=True)
        if not prep_done and clear:
            clear_prep(task_out_dir)
            clear_bin(task_out_dir)
        if stage == 'prep':
            return prep_done

    if stage == 'all' or stage == 'bin':
        task_out_dir = bntsk.task_out_dir
        b = getattr(bn, (bntsk.binner + 'Binner') if not bntsk.is_refiner else (bntsk.refiner+'Refiner'))()
        bin_done =  b.bin_done(**bntsk.to_dict())
        bin_done &= exists(join(task_out_dir, 'bin_DONE'))
        if not bin_done and report:
            print('bin not DONE', bntsk.task_name, flush=True)
        if not bin_done and clear:
            clear_bin(task_out_dir)
        if stage == 'bin':
            return bin_done

    if stage == 'all' or stage == 'qc':
        task_out_dir = bntsk.task_out_dir
        b = getattr(qc, (bntsk.binner + 'Binner') if not bntsk.is_refiner else (bntsk.refiner+'Refiner'))()
        bin_done =  b.bin_done(**bntsk.to_dict())
        bin_done &= exists(join(task_out_dir, 'bin_DONE'))
        if not bin_done and report:
            print('bin not DONE', bntsk.task_name, flush=True)
        if not bin_done and clear:
            clear_bin(task_out_dir)
        if stage == 'bin':
            return bin_done

    #if stage == 'all' or stage == 'post':
    #    spec = json.loads(bntsk.post_spec)
    #    task_out_dir = spec['task_out_dir']
    #    post_done = all(getattr(bn, b+'Binner').post_done(**spec) for b in spec['binners'])
    #    post_done &= exists(join(task_out_dir, 'post_DONE'))
    #    if report:
    #        if only_not_done and not post_done:
    #            print('post_NOT_DONE', bntsk.task_name, bntsk.description)
    #        elif not only_not_done:
    #            print('post_NOT_DONE' if not post_done else 'post_DONE', bntsk.task_name)
    #    if not post_done and clear:
    #        clear_post(spec['task_out_dir'])
    #    if stage == 'post':
    #        return post_done

    return prep_done and bin_done

def bowtie2(idx, bam, r1, r2, threads=8):
    run(f'bowtie2 -p {threads} -x {idx} -1 {r1} -2 {r2} | samtools view -bS - > {bam}', shell=True)

def bowtie2_build(ref, idx):
    run(f'bowtie2-build --threads 4 {ref} {idx}', shell=True,stdout=DEVNULL,stderr=DEVNULL)

def sort_bam(bam, coordinate=True, tmp_pref='tmp', threads=8):
    sort_coordinate = '' if coordinate else ' -n '
    run(f'samtools sort -@ {threads} {sort_coordinate} {bam} > {bam}.{tmp_pref}',shell=True)
    run(f'mv {bam}.{tmp_pref} {bam}',shell=True)

def index_bam(bam):
    bai = bam.replace('.bam', '.bai')
    run(f'samtools index {bam} {bai}',shell=True)

def compute_idxstats(bam, idxstats):
    run(f'samtools idxstats {bam} > {idxstats}',shell=True)