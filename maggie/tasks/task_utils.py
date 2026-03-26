from os.path import join, exists
import pandas as pd
import hashlib
from itertools import product
import json
from ..utils import parse_binning_mode_datasets, manifest_get, write_done_flag, rm_dir
from . import assembly as ab, binning as bn, quality_control as qc

jd = lambda x: json.dumps(x)

def make_manifest(
    assemblers, binners, modes, wrappers, qctools, simulation_dir, 
    assembly_cache, bintask_dir, project_base
):
    bintask_counter = set()
    assembly_counter = set()
    tasks = list()
    for (ab, abo), (bn, bno), mode in product([assemblers, binners, modes]):

        # abo and bno are respectively the assembly options and binning options. 
        # For each option copbination, we make a separate directory
        assembly_counter.add((ab, abo))
        bintask_counter.add((ab, abo, bn, bno, mode))
        spec = {
            'asm_dir': join(assembly_cache, f'assembly_task_{len(assembly_counter)}'),
            'simulation_dir': simulation_dir,
            'task_out_dir': join(bintask_dir, f'bin_task_{len(bintask_counter)}'),
            'assembler': ab, 'assembler_options': abo, 'binner': bn, 'binner_options': bno,
            'binning-mode': mode, 'qctools': qctools,
            'is_wrapper': False
        }
        datasets = parse_binning_mode_datasets(join(project_base, f'{mode}_datasets.csv'))
        for trgt, df in datasets.groupby('target_sample'):
            spec['target_sample'] = trgt
            spec['samples'] = df.dataset.to_list()
            spec['task_out_dir'] = join(spec['task_out_dir'], trgt)
            spec['asm_dir'] = join(spec['asm_dir'], trgt)
            tasks.append(spec)

    for (wrapper, wpo, pipelines) in wrappers:
        bintask_counter.add((wrapper, wpo, pipelines))
        spec = {
            'task_dir': join(bintask_dir, f'bin_task_{len(bintask_counter)}'), 
            'qctools': qctools, 'simulation_dir': simulation_dir, 'is_wrapper': True

        }
        datasets = parse_binning_mode_datasets(join(project_base, f'{mode}_datasets.csv'))
        for trgt, df in datasets.groupby('target_sample'):
            spec['target_sample'] = trgt
            spec['samples'] = df.dataset.to_list()
            spec['task_out_dir'] = join(spec['task_out_dir'], trgt)
            
    for spec in tasks:
        spec['task_name'] = hashlib.sha256(jd(spec.encode()).hexdigest())
    tasks = pd.DataFrame(spec)

def run_assemblies(manifest, threads=8):

    # multiple bintasks can use the sample assembly so first,
    # gather all the unique assembly tasks from the manifest
    asm_tasks = manifest[['target_sample', 'simulation_dir', 'asm_dir', 'assembler', 'assembler_options']].drop_duplicates()
    print('Total assembly tasks: ', len(asm_tasks))
    # run tasks
    for i in range(len(asm_tasks)):
        t = asm_tasks.loc[i,:]
        r1 = join(t.simulation_dir, f'{t.target_sample}_1.fastq.gz')
        r2 = join(t.simulation_dir, f'{t.target_sample}_2.fastq.gz')
        assembler = getattr(ab, f'{t.assembler}Assembler')
        assembler.run_assembly(r1, r2, t.asm_dir, threads, t.assembler_options)
        assembler.clean_up()
        assert assembler.assembly_done()

def run_bin_tasks(manifest, run_only=None, force_prep=False, force_bin=False, threads=8):
    # run prep
    bin_tasks = manifest[~manifest.is_wrapper]
    if run_only is None or run_only == 'prep':
        for i in len(range(bin_tasks)):
            t = bin_tasks.loc[i,:]
            if force_prep or not task_done(t, 'prep', clear=True):
                run_prep(t, threads)
    # run bin
    if run_only is None or run_only == 'bin':
        for i in len(range(manifest)):
            t = manifest.loc[i,:]
            if force_bin or not task_done(t, 'bin', clear=True):
                run_bin(t, threads)

def run_prep(spec, threads):
    print('run_prep:', spec.task_name, flush=True)
    b = getattr(bn, spec.binner + 'Binner')
    b.run_prep(**(spec.to_dict() + {'threads':threads}))
    write_done_flag(spec.task_out_dir, name='prep')

def run_bin(spec, threads):
    print('run_bin:', spec.task_name,flush=True)
    b = getattr(bn, spec.binner + 'Binner')
    b.run_binning(**(spec.to_dict()+ {'threads':threads}))
    write_done_flag(spec.task_out_dir, name='bin')

def clear_prep(task_out_dir):
    rm_dir(join(task_out_dir, 'input'))
    rm_dir(join(task_out_dir, 'prep_DONE'))

def clear_bin(task_out_dir, binner):
    rm_dir(join(task_out_dir, 'output', f'{binner}_bins'))
    rm_dir(join(task_out_dir, 'bin_DONE'))

def task_done(bntsk, stage='all', report=False, only_not_done=False, clear=False):
    if stage == 'all' or stage == 'prep':
        spec = json.loads(bntsk.prep_spec)
        task_out_dir = spec['task_out_dir']
        prep_done = all(getattr(bn, b+'Binner').prep_done(**spec) for b in spec['binners'])
        prep_done &= exists(join(task_out_dir, 'prep_DONE'))
        if report:
            if only_not_done and not prep_done:
                print('prep_NOT_DONE', bntsk.task_name, bntsk.description)
            elif not only_not_done:
                print('prep_NOT_DONE' if not prep_done else 'prep_DONE', bntsk.task_name, bntsk.description)
        if not prep_done and clear:
            clear_prep(spec['task_out_dir'])
            for b in spec['binners']:
                b = getattr(bn, b+'Binner')
                clear_bin(spec['task_out_dir'], b.name)
            clear_post(spec['task_out_dir'])
        if stage == 'prep':
            return prep_done

    if stage == 'all' or stage == 'bin':
        spec = json.loads(bntsk.bin_spec)
        task_out_dir = spec['task_out_dir']
        bin_done = all(getattr(bn, b+'Binner').bin_done(**spec) for b in spec['binners'])
        bin_done &= exists(join(task_out_dir, 'bin_DONE'))
        if report:
            if only_not_done and not bin_done:
                print('bin_NOT_DONE', bntsk.task_name, bntsk.description)
            elif not only_not_done:
                print('bin_NOT_DONE' if not bin_done else 'bin_DONE', bntsk.task_name, bntsk.description)
        if not bin_done and clear:
            for b in spec['binners']:
                b = getattr(bn, b+'Binner')
                clear_bin(spec['task_out_dir'], b.name)
            clear_post(spec['task_out_dir'])
        if stage == 'bin':
            return bin_done

    if stage == 'all' or stage == 'post':
        spec = json.loads(bntsk.post_spec)
        task_out_dir = spec['task_out_dir']
        post_done = all(getattr(bn, b+'Binner').post_done(**spec) for b in spec['binners'])
        post_done &= exists(join(task_out_dir, 'post_DONE'))
        if report:
            if only_not_done and not post_done:
                print('post_NOT_DONE', bntsk.task_name, bntsk.description)
            elif not only_not_done:
                print('post_NOT_DONE' if not post_done else 'post_DONE', bntsk.task_name)
        if not post_done and clear:
            clear_post(spec['task_out_dir'])
        if stage == 'post':
            return post_done

    return prep_done and bin_done and post_done