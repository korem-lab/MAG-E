from os.path import join, exists, basename
import pandas as pd
import hashlib
from functools import reduce
from itertools import product
from collections import defaultdict
import json
from ..utils import rm_dir
from . import assembly as ab, binning as bn, quality_control as qc, coverage as cv 

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

def run_assembly(task, threads=8, force=False, check=False):
    r1 = join(task.simulation_dir, f'{task.target}_R1.fastq.gz')
    r2 = join(task.simulation_dir, f'{task.target}_R2.fastq.gz')
    options = '' if task.assembler_options == 'default' else task.assembler_options
    assembler = getattr(ab, f'{task.assembler}Assembler')()
    kwargs = task.to_dict() | {'threads':threads, 'r1':r1, 'r2':r2,'options':options}
    if check:
        return assembler.main_done(**kwargs)
    if force or not assembler.main_done(**kwargs):
        rm_dir(task.assembler_dir, remake=True)
        assembler.run_main(**kwargs)

def run_coverage(task, stage, cov_sample, threads=8, force=False, check=False):
    if cov_sample is not None:
        assert cov_sample in task.samples
    options = '' if task.coverage_options == 'default' else task.coverage_options
    coverage = getattr(cv, f'{task.coverage}Coverage')()
    kwargs = {'cov_sample':cov_sample, 'threads':threads, 'options':options}
    kwargs = task.to_dict() | kwargs
    if check:
        ckp = coverage.prep_done(**kwargs)
        ckm = coverage.main_done(**kwargs)
        return ckp if stage == 'prep' else ckm if stage == 'main' else (ckp and ckm)
    if stage == 'prep' or stage == 'all':
        if force or not coverage.prep_done(**kwargs):
            rm_dir(task.coverage_dir, remake=True)
            coverage.run_prep(**kwargs)
    if stage == 'main' or stage == 'all':
        if force or not coverage.main_done(**kwargs):
            coverage.run_main(**kwargs)

def run_binning(task, stage, threads=8, force=False, check=False):
    binner = getattr(bn, f'{task.binner}Binner')()
    options = '' if task.binner_options == 'default' else task.binner_options
    kwargs = task.to_dict() | {'threads':threads, 'options':options}
    if check:
        ckp = binner.prep_done(**kwargs)
        ckm = binner.main_done(**kwargs)
        return ckp if stage == 'prep' else ckm if stage == 'main' else (ckp and ckm)
    if stage == 'prep' or stage == 'all':
        if force or not binner.prep_done(**kwargs):
            rm_dir(join(task.binner_dir, 'input'),remake=True)
            binner.run_prep(**kwargs)
    if stage == 'main' or stage == 'all':
        if force or not binner.main_done(**kwargs):
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
#            if not q.has_bins(t.binner_dir):
#                continue
#            if not q.done(t.binner_dir):
#                q.run(**(t.to_dict() | {'threads':threads}))
#            tables.append(q.to_qctable(**(t.to_dict() | {'threads':threads})))
#        if tables:
#            qc_table = reduce(lambda left,right: pd.merge(left,right,on='bin'),tables)
#            qc_table['task_hash'] = t.task_hash
#            qc_table.to_csv(join(t.binner_dir, 'output/qc_table.csv'), index=None)

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
#                    ppln_bin_out.append(bt.binner_dir.item())
#                run_bin(t, threads, refine=True, ppln_bin_out=ppln_bin_out)
