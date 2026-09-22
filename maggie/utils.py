import pandas as pd
import json
import gzip
import ast
from pathlib import Path
import os
from os.path import join
import glob
import shutil
import sys
from subprocess import run as _run

def make_bash_template(name, time, mem, thread):
    return f"""#!/bin/bash 
#SBATCH --job-name={name}
#SBATCH --time={time}
#SBATCH --mem={mem}
#SBATCH --account=pmg
#SBATCH --cpus-per-task={thread}
    """

def get_temp_local(subdir='ic2465'):
    dirname = '/local' if os.path.exists('/local') else '/pmglocal'
    dirname = join(dirname, subdir)
    os.makedirs(dirname, exist_ok=True)
    return dirname

def printit(func):
    def wrapper(*args, **kwargs):
        print(args[0] if args else kwargs.get('cmd'), flush=True)
        return None
    return wrapper
def run(cmd):
    _run(cmd, shell=True)

def add_cmd(cmd, script):
    return script + f"\n{cmd}\n"

def move(src, dest):
    shutil.move(src, dest)

def soft_link(src, dst):
    if not os.path.exists(dst):
        run(f'ln -f -s {src} {dst}')

def print_and_return(s):
    print(s,flush=True,end='')
    return s

def run_R_script(script, *args):
    dir = Path(__file__).parent
    cmd = f'Rscript {dir}/R/{script}.R {" ".join(args)}'
    print(cmd)
    run(cmd)

def write_done_flag(task_out_dir, name):
    Path(os.path.join(task_out_dir, f'{name}_DONE')).touch()

def get_contig_name(hdr):
    """get contig name from fasta header."""
    c = hdr.strip().split()[0]
    if c.startswith('>'):
        c = c[1:]
    return c

def remove_fasta_ext(s):
    return s.replace('.fasta', '').replace('.fa', '')

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
    assert ['sample', 'count'] == df.columns.to_list()
    df.set_index('sample',inplace=True)
    return df

def parse_contig_properties(file):
    df = pd.read_csv(file)
    for e in ['contig', 'sample']:
        assert e in df.columns
    df['key'] = df.apply(lambda x: f'{df.sample}:{df.contig}', axis=1)
    return df

def parse_binning_mode_datasets(file):
    df = pd.read_csv(file)
    assert ['target', 'dataset'] == df.columns.to_list()
    return df

def parse_quality_control_table(file):
    df = pd.read_csv(file)
    return df

def parse_manifest(file):
    df = pd.read_csv(file)
    df.fillna({'assembler_summary_name':'', 'binner_summary_name':'', 'refiner_summary_name':''}, inplace=True)
    df.samples = df.samples.apply(ast.literal_eval)
    if 'qctools' in df.columns:
        df.qctools = df.qctools.apply(ast.literal_eval)
    if 'pipelines' in df.columns:
        df.pipelines = df.pipelines.apply(lambda x: ast.literal_eval(x) if type(x) == str else x)
    return df

def flatten(elem, newl):
    if type(elem) != list:
        newl.append(elem)
    else:
        for subl in elem:
            flatten(subl, newl)

def decompress(genomes, exe='unpigz', script=None):
    cmd = f'{exe} -f {genomes}'
    if script:
        return add_cmd(cmd, script)
    run(cmd)

def compress(genomes, exe='pigz', script=None):
    cmd = f'{exe} -f {genomes}'
    if script:
        return add_cmd(cmd, script)
    run(cmd)

def manifest_get(field, x):
    json.loads(x)[field]

def rm_file(file):
    if os.path.exists(file):
        os.remove(file)
def rm_dir(dir, remake=False):
    try:
        if os.path.exists(dir) and os.path.isdir(dir):
            shutil.rmtree(dir)
            if remake:
                os.makedirs(dir)
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

def not_empty(file):
    return os.path.exists(file) and os.path.getsize(file)>0

def n50(file, gz=False):
    fin = open(file) if not gz else gzip.open(file)
    lens = [len(fa.seq) for fa in parse_fasta(fin, gz)]
    lens = sorted(lens, reverse=True)
    half_sum = sum(lens)/2
    s = 0 
    for i in range(len(lens)):
        s += lens[i]
        if s >= half_sum:
            return lens[i]

def tupleize(obj):
    if isinstance(obj, list):
        return tuple(tupleize(item) for item in obj)
    return obj
