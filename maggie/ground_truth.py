import os
from os.path import join
import pandas as pd
import numpy as np
from subprocess import run
import re
import glob
from .utils import manifest_get, parse_maggie_db

def make_blastdb(ref_genomes, blast_db, title):
    run(f'zcat {ref_genomes} | makeblastdb -dbtype nucl -in - -out {blast_db} -title {title}',shell=True)

def blastn(blast_db, contigs, out_file, threads=8):
    cmd = f'blastn -db {blast_db} -outfmt "6 qacc sacc evalue qstart qend qlen sstart send slen pident nident length" -query {contigs} -out {out_file} -num_threads {threads}'
    run(cmd, shell=True)

def get_blast_file(sample,contigs,genomes, outdir, force=False):
    genomes_str = ' '.join(genomes)
    blast_results = os.path.join(outdir, f'{sample}_blast_results.txt')
    if os.path.exists(blast_results) and not force:
        return load_blast_results(blast_results)
    db_name = os.path.join(outdir, f'bdb_{sample}')
    make_blastdb(genomes_str, db_name, sample)
    blastn(db_name, contigs, blast_results, threads=4)
    return load_blast_results(blast_results)

def filter_blast_hits(hits, min_contig_len=100, min_pident=99, min_prop=99,max_prop=101):
    # TODO: Rather than filtering out contigs, Add a binary, "modelled" column, which is
    # TODO: True if the contig is to be modelled, and false otherwise
    filter = hits.contig_length >= min_contig_len
    # The contig must be a near identical match to the genome
    filter &= hits.pident >= min_pident
    # The length of the match must be nearly the entire length of the contig
    filter &= (min_prop <= hits.aln_length / hits.contig_length * 100) & (hits.aln_length / hits.contig_length * 100 <= max_prop)
    return hits.loc[filter, :]

def load_blast_results(file):
    hits = pd.read_csv(
        file,
        delimiter='\t',
        header=None
    )
    hits.columns = [
            'contig', 'ref', 'evalue', 'contig_start', 'contig_end',
            'contig_length', 'ref_start', 'ref_end', 'ref_length',
            'pident', 'nident', 'aln_length'
        ]
    hits['genome'] = hits.ref.apply(lambda x: x.split('_')[0])
    hits['ref'] = hits.ref.apply(lambda x: x.replace('.fa',''))
    return hits

def _construct_ground_truth(sample, straindb, hits):
    # Initialize ground truth matrix
    isolate_genomes = set(straindb.genome[straindb.GenomeType == 'Isolate'])
    hits['sample'] = sample
    hits.drop(['evalue', 'contig_start', 'contig_end'],axis=1, inplace=True)
    hits['key']  = hits.contig.apply(lambda x : f'{sample}-{x}')
    hits['is_isolate'] = hits.genome.isin(isolate_genomes)
    hits['genome_length'] = hits.genome.map(straindb[['genome', 'Length']].set_index('genome').Length)
    return hits

def retreive_contig_names(hits):
    return sorted(hits.contig.unique())

def construct_ground_truth(
        manifest, maggie_db_table, min_contig_len, min_pident, min_prop, max_prop
    ):
    asm_tasks = set()

    # extract assembly tasks
    asm_tasks = manifest[~manifest.is_refiner][['target_sample', 'simulation_dir', 'assembly_dir']]
    asm_tasks.drop_duplicates(inplace=True)
    
    for i in range(len(asm_tasks)):
        t = asm_tasks.iloc[i,:]
        genomes = glob.glob(join(t.simulation_dir, f'iss_{t.target_sample}_genomes/*.fasta.gz'))
        contigs = join(t.assembly_dir, f'{t.target_sample}.fasta')
        hits = get_blast_file(t.target_sample, contigs, genomes, t.assembly_dir)
        hits = filter_blast_hits(hits, min_contig_len, min_pident, min_prop, max_prop)
        db_table = parse_maggie_db(maggie_db_table)
        gt = _construct_ground_truth(t.target_sample, db_table, hits)
        gt.to_csv(os.path.join(t.assembly_dir, f'{t.target_sample}_gt_table.csv'), index=None)