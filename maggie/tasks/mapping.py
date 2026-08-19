from subprocess import run, DEVNULL
from glob import glob
#    contigs = join(task.asm_dir, f'{task.target_sample}.fasta')
#    bam = contigs.replace('.fasta', f'_{sample}.bam')
#    bowtie2(contigs.replace('.fasta',''),bam, r1, r2, threads)
#    sort_bam(bam,threads=threads)
#    index_bam(bam)

from abc import ABC, abstractmethod
import gzip
from os.path import exists, join, dirname, normpath
from os import rename
import re
from subprocess import run
from ..utils import parse_fasta, rm_dir, rm_file, not_empty

def sort_bam(bam, coordinate=True, tmp_pref='tmp', threads=8):
    sort_coordinate = '' if coordinate else ' -n '
    run(f'samtools sort -@ {threads} {sort_coordinate} {bam} > {bam}.{tmp_pref}',shell=True)
    run(f'mv {bam}.{tmp_pref} {bam}',shell=True)

def index_bam(bam):
    bai = bam.replace('.bam', '.bai')
    run(f'samtools index {bam} {bai}',shell=True)

class Mapper(ABC):
    name: str
    execs: str

    @abstractmethod
    def run_map(self, r1, r2, map_dir, target_sample, map_sample, threads, options, **kwargs):
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
    def map_done(self, map_dir, taregt_sample, map_sample, **kwargs):
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

    def run_map(self, r1, r2, map_dir, target_sample, map_sample, threads, options, **kwargs):
        old_files = glob(f'{map_dir}/*.sam') + glob(f'{map_dir}/*.bam')
        for f in old_files:
            rm_file(f)
        bam = join(map_dir, f'{target_sample}_{map_sample}.bam')
        cmd = f'bowtie2 -p {threads} {options} -x {map_dir}/idx -1 {r1} -2 {r2} | samtools view -bS - > {bam}'
        run(cmd, shell=True)
        sort_bam(bam,threads)
        index_bam(bam)

    def run_prep(self, map_dir, asm_dir, threads, **kwargs):
        rm_dir(asm_dir, remake=True)
        cmd = f'bowtie2-build --threads {threads} {asm_dir}/contigs.fasta {map_dir}/idx'
        run(cmd, shell=True,check=True)
    
    def clean_up(self, map_dir, **kwargs):
        # shouldn't be much to clean with bowtie2, but check for any sams
        files = glob(f'{map_dir}/*.sam')
        for f in files:
            rm_file(f)

    def map_done(self, map_dir, target_sample, map_sample, **kwargs):
        bam = join(map_dir, f'{target_sample}_{map_sample}.bam')
        return not_empty(bam)

    def prep_done(self, map_dir, **kwargs):
        files = glob(f'{map_dir}/*.bt2') + glob(f'{map_dir}/*.bt2l')
        return all(not_empty(f) for f in files)