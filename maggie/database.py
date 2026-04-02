import os
import pandas as pd
from os.path import exists, join
from subprocess import run

def check_genome_files_exist(genomes, directory):
    for g in genomes:
        assert exists(join(directory, f'{g}.fasta.gz'))

def run_strain_clustering(genomes, genome_dir, drep_dir, threads, ani):
    """
    Cluster genomes at the strain level. Currently this is done with dRep.
    """
    # construct the database
    os.makedirs(join(drep_dir), exist_ok=True)
    for sprp, spdf in genomes.groupby('SpeciesRepr'):
         # avoid clustering singlets
         if len(spdf) == 1:
             continue
         gnm_str = ' '.join([f'{genome_dir}/{e}.fasta.gz' for e in spdf.genome])
         cmd = f'dRep dereplicate {drep_dir}/dRep_{sprp} -g {gnm_str} -p {threads} --S_ani {ani} --ignoreGenomeQuality'
         run(cmd, shell=True)

def get_dRep_cluster_info(grep, drep_dir):
    """
    Get the dRep cluster files (Cdb) which shows the strain-level cluster information.
    Construct StrainCID which encodes the strain cluster assignments.
    Using Wdb (which specifies cluster representatives), construct
    isStrainRepr (true, if the genome is the strain representative).
    """
    try:
        cdbfl = join(drep_dir, f'dRep_{grep}/data_tables/Cdb.csv')
        cdb = pd.read_csv(cdbfl)
        cdb.genome = cdb.genome.apply(lambda x: x.replace('.fasta.gz', ''))
        cdb.rename({'secondary_cluster': 'StrainCID'}, axis=1, inplace=True)
        wdbfl = join(drep_dir, f'dRep_{grep}/data_tables/Wdb.csv')
        wdb = pd.read_csv(wdbfl)
        wdb.genome = wdb.genome.apply(lambda x: x.replace('.fasta.gz', ''))
        cdb['isStrainRepr'] = cdb.genome.isin(wdb.genome)
        return cdb
    except FileNotFoundError:
        print('missing', grep)
        return None

def build_database_table(database, genomes_dir, drep_dir):

    # get the genome representatives of species clusters >= 2
    greprs = database.SpeciesRepr.value_counts()[database.SpeciesRepr.value_counts() > 1].index.to_series()

    # get the strain cluster info
    strain_cluster_info = greprs.apply(lambda x: get_dRep_cluster_info(x, drep_dir)).to_list()
    strain_cluster_info = [e for e in strain_cluster_info if e is not None]
    strain_cluster_info = pd.concat(strain_cluster_info)

    # Add strain info to database
    database['SpeciesCID'] = database.SpeciesRepr
    database = pd.merge(database, strain_cluster_info, left_on='genome', right_on='genome', how='outer')
    database['isSpeciesRepr'] = database.genome.isin(database.SpeciesRepr)
    # NA if the species was a singlet, so there was no strain clustering performed. 
    # In this case, the singlet is both the species and strain representatitive, so 
    # set isStrainRepr to true, and set the StrainCID to '1_1'
    database['isStrainRepr'] = database['isStrainRepr'].fillna(True)
    database['StrainCID'] = database['StrainCID'].fillna('1_1')
    # make a StrainCID globally unique by combining the species and current strain CIDs 
    database['StrainCID'] = database.apply(lambda x: f'{x.SpeciesCID}:{x.StrainCID}',axis=1)

    # calculate the strain representatives
    database['StrainRepr'] = database.groupby('StrainCID', group_keys=False).apply(lambda x: pd.Series(x.genome[x.isStrainRepr].values.repeat(len(x)), index=x.index))
    database['FileLocation'] = database.genome.apply(lambda x: f'{genomes_dir}/{x}.fasta.gz')

    return database

def construct_sylphdb(genome_list, maggie_db_dir, db_prefix, c=200, t=8, force=False):
    os.makedirs(maggie_db_dir,exist_ok=True)
    if not force and os.path.exists(f'{maggie_db_dir}/{db_prefix}.syldb'):
        return
    with open(os.path.join(maggie_db_dir, f'{db_prefix}.input_genome_list.txt'), 'w') as f:
        f.write('\n'.join(genome_list))
    run(f'sylph sketch -l {maggie_db_dir}/{db_prefix}.input_genome_list.txt -o {maggie_db_dir}/{db_prefix} -t {t} -c {c}', shell=True)