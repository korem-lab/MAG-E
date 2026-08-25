from subprocess import run, DEVNULL
from glob import glob

from abc import ABC, abstractmethod
import gzip
from os.path import exists, join, dirname, normpath
from os import rename
import re
from ..utils import parse_fasta, rm_dir, rm_file, not_empty, run

def sort_bam(bam, coordinate=True, tmp_pref='tmp', threads=8):
    sort_coordinate = '' if coordinate else ' -n '
    run(f'samtools sort -@ {threads} {sort_coordinate} {bam} > {bam}.{tmp_pref}')
    run(f'mv {bam}.{tmp_pref} {bam}')

def index_bam(bam):
    bai = bam.replace('.bam', '.bai')
    run(f'samtools index {bam} {bai}')

class Mapper(ABC):
    name: str
    execs: str

    @abstractmethod
    def run_main(self, r1, r2, map_dir, target, map_sample, threads, options, **kwargs):
        """
        Runs the core read mapping routine.
        """
        ...

    @abstractmethod
    def run_prep(self, asm_dir, map_dir, threads, **kwargs):
        """
        Runs commands in preparation for mapping (e.g indexing, if needed).
        """
        ...

    @abstractmethod
    def clean_up(self, map_dir, **kwargs):
        """
        Removes unneeded files.
        """
        ...
    @abstractmethod 
    def main_done(self, map_dir, target_sample, map_sample, **kwargs):
        """
        Checks for completed mapping.
        """

    @abstractmethod 
    def prep_done(self, map_dir, **kwargs):
        """
        Checks for completed prep.
        """
    

class bowtie2Mapper(Mapper):
    name = 'MEGAHIT'
    exec = 'megahit'

    def run_main(self, r1, r2, map_dir, target, map_sample, threads, options, **kwargs):
        old_files = glob(f'{map_dir}/*.sam') + glob(f'{map_dir}/*.bam')
        for f in old_files:
            rm_file(f)
        bam = join(map_dir, f'{target}_{map_sample}.bam')
        cmd = f'bowtie2 -p {threads} {options} -x {map_dir}/idx -1 {r1} -2 {r2} | samtools view -bS - > {bam}'
        run(cmd)
        sort_bam(bam,threads)
        index_bam(bam)

    def run_prep(self, map_dir, asm_dir, threads, **kwargs):
        rm_dir(asm_dir, remake=True)
        cmd = f'bowtie2-build --threads {threads} {asm_dir}/contigs.fasta {map_dir}/idx'
        run(cmd)
    
    def clean_up(self, map_dir, **kwargs):
        pass

    def main_done(self, map_dir, target, map_sample, **kwargs):
        bam = join(map_dir, f'{target}_{map_sample}.bam')
        return not_empty(bam)

    def prep_done(self, map_dir, **kwargs):
        files = glob(f'{map_dir}/*.bt2') + glob(f'{map_dir}/*.bt2l')
        return all(not_empty(f) for f in files) and files