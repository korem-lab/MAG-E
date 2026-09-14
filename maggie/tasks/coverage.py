from subprocess import run, DEVNULL
from glob import glob
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

class Bowtie2JGI(Coverage):

    def run_main(self, cov_dir, abund_dir, target, samples, **kwargs):
        # MAKE DIFFERENT FORMATS 
        pass

    def run_prep(self, cov_dir, target, samples, abund_dir, **kwargs):
        bams = [join(cov_dir, f'{target}_{s}.bam') for s in samples]
        jgi_summarize(bams, join(abund_dir, f'coverage.tsv'))

    def main_done(self, abund_dir, **kwargs):
        # FORMAT CHECK
        pass

    def prep_done(self, abund_dir, **kwargs):
        return not_empty(join(abund_dir, f'coverage.tsv'))

class BWAMEM2JGI(Coverage):

    def run_main(self, cov_dir, abund_dir, target, samples, **kwargs):
        # MAKE DIFFERENT FORMATS 
        pass

    def run_prep(self, cov_dir, target, samples, abund_dir, **kwargs):
        bams = [join(cov_dir, f'{target}_{s}.bam') for s in samples]
        jgi_summarize(bams, join(abund_dir, f'coverage.tsv'))

    def main_done(self, abund_dir, **kwargs):
        # FORMAT CHECK
        pass

    def prep_done(self, abund_dir, **kwargs):
        return not_empty(join(abund_dir, f'coverage.tsv'))

class Fairy(Coverage):
    name: 'fairy'
    execs: 'fairy'

    def run_main(self, **kwargs):
        # MAKE ALL FORMATS AVAILABLE HERE
        pass

    def run_prep(self, asm_dir, simulation_dir, abund_dir, samples, threads, **kwargs):
        # make sketches
        for s in samples:
            s = join(simulation_dir, s)
            cmd = f'{self.name} sketch 1 {s}_R1.fastq.gz -2 {s}_R2.fastq.gz -d {abund_dir}'
            run(cmd)

        # calculate coverage matrix
        contigs = join(asm_dir, 'contigs.fasta')
        cmd = f'{self.name} coverage {abund_dir}/*.bcsp {contigs} -t {threads} -o {abund_dir}/coverage_jgifmt.tsv'
        run(cmd)
        cmd = f'{self.name} coverage --maxbin-format {abund_dir}/*.bcsp {contigs} -t {threads} -o {abund_dir}/coverage_mxbfmt.tsv'
        run(cmd)

        # calculate individual coverage files for SemiBin2
        cmd = f'SemiBin2 split_contigs -i {contigs} -o {abund_dir}/semibin2_splits'
        run(cmd)
        for s in samples:
            s = join(simulation_dir, s)
            cmd = f'{self.name} coverage {abund_dir}/semibin2_splits/split_contigs.fna.gz {abund_dir}/{s}.paired.bcsp' \
            f'--aemb-format -o {abund_dir}/{s}_coverage_aembfmt.tsv'
            run(cmd)

        #cleanup
        rm_dir(f'{abund_dir}/semibin2_splits')
        for s in samples:
            rm_file(f'{abund_dir}/{s}.bcsp')


    def main_done(self, abund_dir, samples, **kwargs):
        # CHECK ALL FORMATS AVAILABLE
        pass

    def prep_done(self, abund_dir, samples, **kwargs):
        return not_empty(join(abund_dir, 'coverage.tsv')) and \
            all(not_empty(join(abund_dir, f'{s}_coverage.tsv')) for s in samples)

class AEMB(Coverage):
    name: 'strobealign'
    execs: 'strobealign'

    def run_main(self, simulation_dir, asm_dir, abund_dir, target, cov_sample, threads, **kwargs):
        # FORMAT FOR DIFFERENT BINNERS
        pass 

    def run_prep(self, simulation_dir, asm_dir, abund_dir, cov_sample, threads, **kwargs):
        contigs = join(asm_dir, 'contigs.fasta')
        r1 = join(simulation_dir, f'{cov_sample}_R1.fastq.gz')
        r2 = join(simulation_dir, f'{cov_sample}_R2.fastq.gz')
        cmd = f'{self.execs} -t {threads} --aemb {contigs} {r1} {r2} > {abund_dir}/{cov_sample}_coverage.tsv'
        run(cmd)

    def main_done(self, cov_dir, target_sample, cov_sample, **kwargs):
        # CHECK ALL FORMATS AVAILABLE
        pass

    def prep_done(self, abund_dir, cov_sample, **kwargs):
        return not_empty(join(abund_dir, f'{cov_sample}_coverage.tsv'))
