import os
import glob
import pandas as pd
from os.path import join, exists
from subprocess import run
from .utils import parse_read_counts, decompress, compress

def run_InSilicoSeq(spec, simout, read_counts, prefix1, n_reads='auto', threads=8, force=True):

    # Either supply a set number of reads per sample
    # Or supply 'auto' and the number of reads in each simulation will match the number of reads in the real data
    if n_reads == 'auto':
        read_counts = parse_read_counts(read_counts)
    else:
        assert type(n_reads) == int

    # copy genomes to a sample-specific directory 
    def copy_over_genomes(x):
        y = join(genomes_dir, x.split('/')[-1])
        run(f'cp {x} {y}',shell=True)
    spec = pd.read_csv(spec)
    sample = spec.Sample_file.iloc[0]
    genomes_dir = join(simout, f'iss_{sample}_genomes')
    os.makedirs(genomes_dir, exist_ok=True)
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
    prefix2 = prefix1.replace('1', '2')
    if not force and (exists(f'{out_pref}_{prefix1}.fastq.gz') and exists(f'{out_pref}_{prefix2}.fastq.gz')):
        print(f'Skipping {sample}. Simulated reads already exist. Delete if you want to re-simulate.', flush=True)
        return
    decompress(' '.join(e+'.gz' for e in list(spec.Genome_file)))
    genomes_str = ' '.join(list(spec.Genome_file))
    _n_reads  = 2*(n_reads if (n_reads != "auto") else read_counts.loc[sample, 'count'])
    run(f'iss generate -p {threads} --draft {genomes_str} --abundance_file {abundance_file} -o {out_pref} --model HiSeq --n_reads {int(_n_reads)} --seed {SEED}',shell=True)
    reads = glob.glob(f'{out_pref}*.fastq')
    print('Compressing reads and genomes...', flush=True)
    compress(genomes_str)
    compress(' '.join(reads))
