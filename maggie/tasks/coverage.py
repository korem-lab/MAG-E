from subprocess import run, DEVNULL
from glob import glob
import pandas as pd
from abc import ABC, abstractmethod
from os.path import exists, splitext, join
from ..utils import rm_dir, rm_file, not_empty, run
from . import assembly as ab

def sort_bam(bam, coordinate=True, tmp_pref='tmp', threads=8):
    sort_coordinate = '' if coordinate else ' -n '
    run(f'samtools sort -@ {threads} {sort_coordinate} {bam} > {bam}.{tmp_pref}')
    run(f'mv {bam}.{tmp_pref} {bam}')

def index_bam(bam):
    bai = bam.replace('.bam', '.bai')
    run(f'samtools index {bam} {bai}')


def jgi_summarize(bams, filename):
    log = splitext(filename)[0] + '.log'
    run(f'jgi_summarize_bam_contig_depths --outputDepth {filename} {bams} 2> {log}')

def jgi_to_maxbin2(dir, coverage):
    df = pd.read_csv(coverage, delimiter='\t')
    for s in df.columns[3::2]:
        df.to_csv[['contigName',s]].to_csv(
            join(dir, f'{s}_coverage.tsv'),
            index=None, header=None, sep='\t'
        )

class Coverage(ABC):
    name: str
    execs: str

    @abstractmethod
    def run_main(self, r1, r2, cov_dir, target, cov_sample, threads, options, **kwargs):
        """
        Runs the core read coverage routine.
        """
        ...

    @abstractmethod
    def run_prep(self, asm_dir, cov_dir, threads, **kwargs):
        """
        Runs commands in preparation for coverage (e.g indexing, if needed).
        """
        ...

    @abstractmethod
    def clean_up(self, cov_dir, **kwargs):
        """
        Removes unneeded files.
        """
        ...

    @abstractmethod 
    def main_done(self, cov_dir, target_sample, cov_sample, **kwargs):
        """
        Checks for completed mapping.
        """
        ...

    @abstractmethod 
    def prep_done(self, cov_dir, **kwargs):
        """
        Checks for completed prep.
        """
        ...

class bowtie2Coverage(Coverage):
    name = 'MEGAHIT'
    exec = 'megahit'

    def run_main(self, r1, r2, cov_dir, target, cov_sample, threads, options, **kwargs):
        old_files = [
            f'{cov_dir}/{target}_{cov_sample}.sam',
            f'{cov_dir}/{target}_{cov_sample}.bam',
        ]
        for f in old_files:
            rm_file(f)
        bam = join(cov_dir, f'{target}_{cov_sample}.bam')
        cmd = f'bowtie2 -p {threads} {options} -x {cov_dir}/idx -1 {r1} -2 {r2} | samtools view -bS - > {bam}'
        run(cmd)
        sort_bam(bam,threads)
        index_bam(bam)

    def run_prep(self, cov_dir, asm_dir, threads, **kwargs):
        rm_dir(cov_dir, remake=True)
        cmd = f'bowtie2-build --threads {threads} {asm_dir}/contigs.fasta {cov_dir}/idx'
        run(cmd)
    
    def clean_up(self, cov_dir, **kwargs):
        pass

    def main_done(self, cov_dir, target, cov_sample, **kwargs):
        bam = join(cov_dir, f'{target}_{cov_sample}.bam')
        return not_empty(bam)

    def prep_done(self, cov_dir, **kwargs):
        # check the assembly is present
        files = glob(f'{cov_dir}/*.bt2') + glob(f'{cov_dir}/*.bt2l')
        return all(not_empty(f) for f in files) and files

class bowtie2JGI(Coverage):

    def run_main(self, cov_dir, **kwargs):
        # get to the coverage file 
        bams = glob(join(cov_dir, '*.bam'))
        jgi_summarize(bams, join(cov_dir, 'coverage.tsv'))

        # make MaxBin2 format
        jgi_to_maxbin2(cov_dir, join(cov_dir, 'coverage.tsv'))


    def run_prep(self, cov_dir, simulation_dir, target, samples, asm_dir, threads, options, **kwargs):
        # remove old files
        rm_dir(cov_dir, remake=True)
        for f in glob.glob(join(cov_dir, '*.bam')):
            rm_file(f)

        # index
        cmd = f'bowtie2-build --threads {threads} {asm_dir}/contigs.fasta {cov_dir}/idx'
        run(cmd)

        # map
        for s in samples:
            r1 = join(simulation_dir, f'{s}_R1.fastq.gz')
            r2 = join(simulation_dir, f'{s}_R2.fastq.gz')
            bam = join(cov_dir, f'{target}_{s}.bam')
            cmd = f'bowtie2 -p {threads} {options} -x {cov_dir}/idx -1 {r1} -2 {r2} | samtools view -bS - > {bam}'
            run(cmd)
            sort_bam(bam, threads)
            index_bam(bam)

    def main_done(self, cov_dir, **kwargs):
        fls = glob.glob(join(cov_dir, '*.tsv'))
        return all(not_empty(e) for e in fls)

    def prep_done(self, cov_dir, target, samples, **kwargs):
        return all(not_empty(join(cov_dir, f'{target}_{s}.bam')) for s in samples)

class bwamem2JGI(Coverage):

    def run_main(self, cov_dir, **kwargs):
        # get to the coverage file 
        bams = glob(join(cov_dir, '*.bam'))
        jgi_summarize(bams, join(cov_dir, 'coverage.tsv'))

        # make MaxBin2 format
        jgi_to_maxbin2(cov_dir, join(cov_dir, 'coverage.tsv'))


    def run_prep(self, cov_dir, simulation_dir, target, samples, asm_dir, threads, options, **kwargs):
        # remove old files
        rm_dir(cov_dir, remake=True)
        for f in glob.glob(join(cov_dir, '*.bam')):
            rm_file(f)

        # index
        cmd = f'bwa-mem2 index --threads {threads} -p {cov_dir}/idx {asm_dir}/contigs.fasta '
        run(cmd)

        # map
        for s in samples:
            r1 = join(simulation_dir, f'{s}_R1.fastq.gz')
            r2 = join(simulation_dir, f'{s}_R2.fastq.gz')
            bam = join(cov_dir, f'{target}_{s}.bam')
            cmd = f'bwa-mem2 mem -t {threads} {options}  {cov_dir}/idx {r1} {r2} | samtools view -bS - > {bam}'
            run(cmd)
            sort_bam(bam, threads)
            index_bam(bam)

    def main_done(self, cov_dir, **kwargs):
        fls = glob.glob(join(cov_dir, '*.tsv'))
        return all(not_empty(e) for e in fls)

    def prep_done(self, cov_dir, target, samples, **kwargs):
        return all(not_empty(join(cov_dir, f'{target}_{s}.bam')) for s in samples)

class Fairy(Coverage):
    name: 'fairy'
    execs: 'fairy'

    def run_main(self, asm_dir, cov_dir, threads, **kwargs):
        # calculate coverage matrix
        contigs = join(asm_dir, 'contigs.fasta')
        cmd = f'fairy coverage {cov_dir}/*.bcsp {contigs} -t {threads} -o {cov_dir}/coverage.tsv'
        run(cmd)
        jgi_to_maxbin2(cov_dir, join(cov_dir, 'coverage.tsv'))

    def run_prep(self, simulation_dir, cov_dir, samples, threads, **kwargs):
        # make sketches
        for s in samples:
            s = join(simulation_dir, s)
            cmd = f'fairy sketch -t {threads} -1 {s}_R1.fastq.gz -2 {s}_R2.fastq.gz -d {cov_dir}'
            run(cmd)

        ## calculate individual coverage files for SemiBin2
        #cmd = f'SemiBin2 split_contigs -i {contigs} -o {abund_dir}/semibin2_splits'
        #run(cmd)
        #for s in samples:
        #    s = join(simulation_dir, s)
        #    cmd = f'{self.name} coverage {abund_dir}/semibin2_splits/split_contigs.fna.gz {abund_dir}/{s}.paired.bcsp' \
        #    f'--aemb-format -o {abund_dir}/{s}_coverage_aembfmt.tsv'
        #    run(cmd)

        ##cleanup
        #rm_dir(f'{abund_dir}/semibin2_splits')
        #for s in samples:
        #    rm_file(f'{abund_dir}/{s}.bcsp')


    def main_done(self, cov_dir, **kwargs):
        fls = glob.glob(join(cov_dir, '*.tsv'))
        return all(not_empty(e) for e in fls)

    def prep_done(self, cov_dir, samples, **kwargs):
        return all(not_empty(join(cov_dir, f'{s}.bcsp')) for s in samples)

class AEMB(Coverage):
    name: 'strobealign'
    execs: 'strobealign'

    def run_main(self, **kwargs):
        pass 

    def run_prep(self, simulation_dir, asm_dir, cov_dir, samples, threads, **kwargs):
        contigs = join(asm_dir, 'contigs.fasta')
        for s in samples:
            r1 = join(simulation_dir, f'{s}_R1.fastq.gz')
            r2 = join(simulation_dir, f'{s}_R2.fastq.gz')
            cmd = f'{self.execs} -t {threads} --aemb {contigs} {r1} {r2} > {cov_dir}/{s}_coverage.tsv'
            run(cmd)

    def main_done(self, cov_dir, target_sample, cov_sample, **kwargs):
        pass

    def prep_done(self, abund_dir, cov_sample, **kwargs):
        return not_empty(join(abund_dir, f'{cov_sample}_coverage.tsv'))
