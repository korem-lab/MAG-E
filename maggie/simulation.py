import os
import glob
import pandas as pd
from os.path import join, exists,basename
from subprocess import run
import numpy as np
from .utils import parse_read_counts, decompress, compress, make_bash_template, add_cmd

def run_InSilicoSeq(spec, simout, read_counts, n_reads='auto', threads=8, seed=37, force=True, print_script=False):
    np.random.seed(37)
    # Either supply a set number of reads per sample
    # Or supply 'auto' and the number of reads in each simulation will match the number of reads in the real data
    if n_reads == 'auto':
        read_counts = parse_read_counts(read_counts)
    else:
        n_reads = int(n_reads)

    # copy genomes to a sample-specific directory 
    def copy_over_genomes(x):
        y = join(genomes_dir, x.split('/')[-1])
        run(f'cp {x} {y}',shell=True)
    spec = pd.read_csv(spec)
    sample = basename(spec.Sample_file.iloc[0])
    genomes_dir = join(simout, f'iss_{sample}_genomes')
    os.makedirs(genomes_dir, exist_ok=True)
    newloc = '/insomnia001/depts/pmg/KoremLab/Projects/MAG-E_analysis/data/ecosystem_databases/UHGG/genomes'
    spec.Genome_file = spec.Genome_file.apply(lambda x: newloc + '/' + x.split('/')[-1])
    spec.Genome_file.apply(copy_over_genomes)
    spec.Genome_file = spec.Genome_file.apply(
        lambda x: join(genomes_dir, os.path.basename(x))
    )

    # make abundance spec
    abundance_file = join(simout, f'{sample}_abundance.txt')
    spec.Genome_file = spec.Genome_file.apply(lambda x: x.replace('.gz', ''))
    spec[['Genome_file','StrainAbund']].to_csv(abundance_file, index=None, header=None, sep='\t')

    # define output prefix for simulation
    out_pref = join(simout, f'{sample}')

    # run iss
    if not force and (exists(f'{out_pref}_R1.fastq.gz') and exists(f'{out_pref}_R2.fastq.gz')):
        print(f'Skipping {sample}. Simulated reads already exist. Delete if you want to re-simulate.', flush=True)
        return
    script = make_bash_template(name=sample, time='12:0:0', mem='32G', thread=8)
    script = decompress(' '.join(e+'.gz' for e in list(spec.Genome_file)), script=script)
    genomes_str = ' '.join(list(spec.Genome_file))
    _n_reads  = 2*(n_reads if (n_reads != "auto") else read_counts.loc[sample, 'count'])
    cmd = f'iss generate -p {threads} --draft {genomes_str} --abundance_file {abundance_file} -o {out_pref} --model HiSeq --n_reads {int(_n_reads)} --seed {seed}'
    if print_script:
        script = add_cmd(cmd, script)
    else:
        run(cmd, shell=True)
    #iss has a bug where temp files are left behind.
    slptm = 300
    cmd = f'sleep {slptm}'

    if print_script:
        script = add_cmd(cmd, script)
        script = compress(genomes_str, script=script)
        script = add_cmd(f'echo DONE > {simout}/{sample}_DONE', script)
        with open(f'{out_pref}_runiss.sh','w') as f:
            f.write(script + '\n')
        return
    print(f'Sleeping {slptm}s to wait for InSilicoSeq to cleanup...', flush=True)
    run(cmd, shell=True)

    tmp_files = glob.glob(f'{out_pref}.iss.tmp*')
    for f in tmp_files:
        os.remove(f)
    reads = glob.glob(f'{out_pref}*.fastq')
    print('Compressing reads and genomes...', flush=True)
    compress(genomes_str)
    compress(' '.join(reads))
