import glob
import pandas as pd
import numpy as np
from pandas.api.types import CategoricalDtype
from os.path import join, basename
from os import makedirs
from subprocess import run

def construct_sylphsp(reads_dir, sylphsp_dir, prefix1, c=200, t=8):
    """
    For all paired-end samples within a directory (must be gziped with .fastq.gz prefix)
    construct a Sylph Sample sketch (.paired.sp).
    """
    makedirs(sylphsp_dir, exist_ok=True)
    prefix2 = prefix1.replace('1', '2')
    r1_reads = sorted(glob.glob(join(reads_dir, f'*_{prefix1}.fastq.gz')))
    r2_reads = sorted(glob.glob(join(reads_dir, f'*_{prefix2}.fastq.gz')))
    sylph_bases = [
        join(sylphsp_dir, e.replace(f'_{prefix1}.fastq.gz', '').split('/')[-1]) for e in r1_reads
    ]
    r1_reads = ' '.join(r1_reads)
    r2_reads = ' '.join(r2_reads)
    bases_str = ' '.join(sylph_bases)
    cmd = f'sylph sketch -1 {r1_reads} -2 {r2_reads} -S {bases_str} -c {c} -t {t} -d {sylphsp_dir}'
    run(cmd, shell=True)
    return [f'{e}.paired.sylsp' for e in sylph_bases]

def sylph_profile(sylsp_dir, maggie_db_dir, sim_dir, db_prefix, t=8):
    run(f'sylph profile {sylsp_dir}/*.sylsp {maggie_db_dir}/{db_prefix}.syldb -t {t} -o {sim_dir}/sylph_profile.tsv',shell=True)
    # make per-sample files for easier loading
    for sample, df in pd.read_csv(f'{sim_dir}/sylph_profile.tsv', sep='\t').groupby('Sample_file'):
        sample_prof = join(sim_dir, f'{basename(sample)}_sylph_profile.tsv')
        df.to_csv(sample_prof, sep='\t', index=None)
    return join(sim_dir, 'sylph_profile.tsv')

def sylph_query(sylsp_dir, maggie_db_dir, sim_dir, db_prefix, t=8):
    run(f'sylph query --minimum-ani 95 {sylsp_dir}/*.sylsp {maggie_db_dir}/{db_prefix}.syldb -t {t} -o {sim_dir}/sylph_query.tsv', shell=True)
    # make per-sample files for easier loading
    for sample, df in pd.read_csv(f'{sim_dir}/sylph_query.tsv', sep='\t').groupby('Sample_file'):
        sample_prof = join(sim_dir, f'{basename(sample)}_sylph_query.tsv')
        df.to_csv(sample_prof, sep='\t', index=None)
    return join(sim_dir, 'sylph_query.tsv')

def select_genomes(df, ani=99.8, ci_lower=99.5, max_strains=3):

    # Select all genomes that have ANI >= 99.8 with either a high confidence adjusted ANI or naive ANI
    fltr = (df.Adjusted_ANI >= ani) & (df['ANI_5-95_percentile'].apply(lambda x: True if 'NA-NA'== x else float(x.split('-')[0]) >= ci_lower))
    df_filt = df[fltr]

    if not df_filt.empty:
        # We define multiple matches of this quality where each match is in a different
        # tertiary cluster as evidence of multiple strains. In this case, pick a random 
        # genome from each match in a different tertiary cluster to replicate strain variability. 
        # Prefer isolates when possible.
        df_filt = df_filt.groupby('StrainCID',group_keys=False).apply(lambda x: x.sort_values(['GenomeType', 'N50', 'Adjusted_ANI'], ascending=False).iloc[0]).reset_index(drop=True)
        if max_strains is not None and len(df_filt) > max_strains:
            df_filt = df_filt.iloc[:max_strains,:]
        # Should multiple strains be present, set their abundances to a log-normal distribtion, summing to the species level abundance
        abn = np.random.lognormal(mean=1, sigma=2, size=len(df_filt))
        abn = abn/abn.sum()
        df_filt['StrainAbund'] = df_filt.SpeciesAbund * abn
        df_filt['StrainAbundUnscaled'] = abn
        return df_filt
    else:
        # if nothing matches at 99.8 with high confidence, take the single genome with the highest ANI
        df = df.sort_values(['GenomeType', 'Adjusted_ANI', 'N50'], ascending=False).iloc[0]
        df['StrainAbund'] = df.SpeciesAbund
        df['StrainAbundUnscaled'] = 1.0
        return pd.DataFrame(df).T

def construct_metagenomic_specification(sim_dir, db_table, prof, query):
    # load query and profile
    prof = pd.read_csv(prof, delimiter='\t')
    query = pd.read_csv(query, delimiter='\t')
    prof['genome'] = prof.Genome_file.apply(lambda x: x.split('/')[-1].replace('.fasta.gz', ''))
    query['genome'] = query.Genome_file.apply(lambda x: x.split('/')[-1].replace('.fasta.gz', ''))

    # add cluster info to the query ANI measurements 
    db_table = pd.read_csv(db_table)
    cat_dtype = CategoricalDtype(categories=['MAG', 'Isolate'], ordered=True)
    db_table.GenomeType = db_table.GenomeType.astype(cat_dtype)
    query = query.groupby('Sample_file',group_keys=False).apply(lambda qry: pd.merge(qry, db_table, on='genome')).reset_index()
    query.drop('index', axis=1, inplace=True)

    # Add the abundance information:
    # - Consider each abundance of each genome in a species cluster to be the same. 
    # - We outerjoin the abundance information in the profile with the cluster info from db_table now merged with the sylph-query.
    genome_abundances = pd.merge(query, prof[['Sample_file', 'genome', 'Sequence_abundance', 'Taxonomic_abundance', 'kmers_reassigned']], on=['Sample_file', 'genome'], how='outer')
    genome_abundances.rename({'Sequence_abundance':'SpeciesAbund'},axis=1,inplace=True)
    if any(genome_abundances.Adjusted_ANI.isna()):
        print('WARNING: following were in the profile but not in the query:')
        print(genome_abundances[genome_abundances.Adjusted_ANI.isna()])
    # Back and forward fill to copy the abundance info from the representative to the entire species cluster 
    genome_abundances= genome_abundances.groupby(['Sample_file', 'SpeciesCID'],group_keys=False).apply(lambda x: x.bfill().ffill()).reset_index()
    genome_abundances.drop('index', axis=1, inplace=True)

    # Sylph profile and query can have discrepancy.
    # Specifically, species-representatives can have ANI > 95 in query, but be absent from the
    # profile, despite the profile also being at ANI > 95. We noticed this typically applies to
    # genomes whose 1) ANI is close to the bound of 95, and 2) where sylph has no confidence in the
    # reported ANI score.
    # In addition to this, it could be that even though the representitive is not in the profile,
    # many of the cluster's genomes have ANI > 95. 
    # In both these cases (representitive in query, but not in profile, whilst cluster's genomes are in query OR representative not in query, also not in profile, whilst cluster genomes are in query)
    # there is no SpeciesAbund to populate to the rest of the cluster. 
    # All such genomes can be dropped - we have no way of quantifying their abundance.
    genome_abundances = genome_abundances[genome_abundances.SpeciesAbund.notna()]

    # We also remove all genomes for which sylph could not compute EffLambda, because coverage was too low.
    # In this case, neither adjusted nor naive ANI will be able to accurately calculate abundance, so drop them from analysis. 
    genome_abundances = genome_abundances[(genome_abundances['Eff_lambda'] != 'LOW')]
    # renormalize after filtering, set to 0-1 scale.
    genome_abundances.SpeciesAbund /= genome_abundances[genome_abundances.isSpeciesRepr].SpeciesAbund.sum()

    # construct metagenomic specifications
    specs = genome_abundances.groupby(['Sample_file', 'SpeciesCID'],group_keys=False).apply(select_genomes).reset_index(drop=True)
    # Set range to 0-1 (rather than 0-100)

    # set paths to absolute
    sample = basename(specs.Sample_file.iloc[0])
    specs.to_csv(join(sim_dir, f'{sample}_metagenome_spec.csv'), index=None)
