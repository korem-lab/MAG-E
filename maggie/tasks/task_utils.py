from os.path import join, exists
import pandas as pd
import hashlib
from itertools import product
import json
from ..utils import parse_binning_mode_datasets, manifest_get
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
            'binning-mode': mode, 'qctools': qctools
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
            'qctools': qctools, 'simulation_dir': simulation_dir

        }
        datasets = parse_binning_mode_datasets(join(project_base, f'{mode}_datasets.csv'))
        for trgt, df in datasets.groupby('target_sample'):
            spec['target_sample'] = trgt
            spec['samples'] = df.dataset.to_list()
            spec['task_out_dir'] = join(spec['task_out_dir'], trgt)
            tasks.append(spec)
    tasks = pd.DataFrame(
        [
            (hashlib.sha256(jd(e.encode()).hexdigest()), jd(e)) for e in tasks
        ], 
        columns=['task_name', 'spec']
    )
    return tasks

def run_assemblies(manifest, threads=8):

    # multiple bintasks can use the sample assembly so first,
    # gather all the unique assembly tasks from the manifest
    asm_tasks = set()
    for i in range(len(manifest)):
        spec = manifest.loc[i, 'spec']
        ts = manifest_get(spec, 'target_sample')
        sd = manifest_get(spec, 'simulation_dir')
        r1 = join(sd, f'{ts}_1.fastq.gz')
        r2 = join(sd, f'{ts}_2.fastq.gz')
        asm_dir = manifest_get(spec, 'asm_dir')
        asm_name = manifest_get(spec, 'assembler')
        asm_options = manifest_get(spec, 'assembler_options')
        asm_tasks.add( 
            (r1, r2, asm_dir, asm_name, asm_options) 
        )
    print('Total assembly tasks: ', len(asm_tasks))
    # run tasks
    for (r1, r2, asm_dir, asm_name, asm_options) in asm_tasks:
        assembler = getattr(ab, f'{asm_name}Assembler')
        assembler.run_assembly(r1, r2, asm_dir, threads, asm_options)
        assembler.clean_up()
        assert assembler.assembly_done()