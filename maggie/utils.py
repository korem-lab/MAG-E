import pandas as pd
import json
from pathlib import Path
import os
import shutil
import sys
from subprocess import run

def write_done_flag(task_out_dir, name):
    Path(os.path.join(task_out_dir, f'{name}_DONE')).touch()

def parse_maggie_db(file):
    df = pd.read_csv(file)
    for e in [
        'StrainCID', 'SpeciesCID', 'isSpeciesRepr', 
        'isStrainRepr', 'StrainRepr', 'FileLocation'
    ]:
        assert e in df.columns
    return df
        
def parse_read_counts(file):
    df = pd.read_csv(file)
    assert ['samples', 'count'] == df.columns
    df.set_index('samples',inplace=True)
    return df

def parse_binning_mode_datasets(file):
    df = pd.read_csv(file)
    assert ['target_sample', 'dataset'] == df.columns
    return df

def parse_manifest(file):
    df = pd.read_csv(file):
    assert ['task_name', 'spec'] == df.columns
    return df

def decompress(genomes, exe='pigz'):
    run(f'{exe} -f {genomes}', shell=True)

def compress(genomes, exe='pigz'):
    run(f'{exe} -f {genomes}', shell=True)

def manifest_get(field, x):
    json.loads(x)[field]


def rm_dir(dir, remake=False):
    try:
        if os.path.exists(dir):
            if os.path.isdir(dir):
                shutil.rmtree(dir)
            elif os.path.isfile(dir):
                os.remove(dir)
    except FileNotFoundError:
        if remake:
            os.makedirs(dir)
    except OSError as e:
        print(f'error deleting {dir}: {e}')
        sys.exit()

class Fasta:
    def __init__(self, hdr, seq, fold=80):
        self.hdr = hdr
        self.seq = seq
        self.fold = fold

    def write(self, fs):
        fs.write(
            f'>{self.hdr}\n' +
            '\n'.join([self.seq[i:i+self.fold] for i in range(0, len(self.seq), self.fold)]) +
            '\n'
        )

def parse_fasta(fs, gzipped=False):
    '''Input: list of contig fasta file names.'''

    decode = lambda x: x if not gzipped else x.decode('utf-8')
    line = decode(next(fs))
    while line[0] == '>':
        # strip header, get rid of '>', and split by ':'
        hdr = line.strip()[1:]
        line = decode(next(fs))  # move to seq
        seq = list()
        while line[0] != '>':
            seq.append(line.strip())
            line = next(fs, False)
            if not line:
                break
            line = decode(line)
        # reached header for next elem
        yield Fasta(hdr, ''.join(seq))
        if not line:
            break