import pandas as pd
from subprocess import run

def parse_read_counts(read_counts):
    df = pd.read_csv(read_counts)
    assert ['samples', 'count'] == df.columns
    df.set_index('samples',inplace=True)
    return df

def decompress(genomes, exe='pigz'):
    run(f'{exe} -f {genomes}', shell=True)

def compress(genomes, exe='pigz'):
    run(f'{exe} -f {genomes}', shell=True)