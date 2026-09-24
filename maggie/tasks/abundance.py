from abc import ABC, abstractmethod
from os.path import exists, splitext, join
import glob
from ..utils import run, not_empty, rm_dir, rm_file

def jgi_summarize(bams, filename):
    log = splitext(filename)[0] + '.log'
    run(f'jgi_summarize_bam_contig_depths --outputDepth {filename} {bams} 2> {log}')

class Coverage(ABC):
    name: str
    execs: str

    @abstractmethod
    def run_main(self, r1, r2, map_dir, target, map_sample, threads, options, **kwargs):
        """
        Runs the core read mapping routine.
        """
        ...

    @abstractmethod
    def run_prep(self, assembly_dir, map_dir, threads, **kwargs):
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
        ...

    @abstractmethod 
    def prep_done(self, map_dir, **kwargs):
        """
        Checks for completed prep.
        """
        ...


class Bowtie2JGI(Coverage):

    def run_main(self, map_dir, abund_dir, target, samples, **kwargs):
        # MAKE DIFFERENT FORMATS 
        pass

    def run_prep(self, map_dir, target, samples, abund_dir, **kwargs):
        bams = [join(map_dir, f'{target}_{s}.bam') for s in samples]
        jgi_summarize(bams, join(abund_dir, f'coverage.tsv'))

    def main_done(self, abund_dir, **kwargs):
        # FORMAT CHECK
        pass

    def prep_done(self, abund_dir, **kwargs):
        return not_empty(join(abund_dir, f'coverage.tsv'))

class BWAMEM2JGI(Coverage):

    def run_main(self, map_dir, abund_dir, target, samples, **kwargs):
        # MAKE DIFFERENT FORMATS 
        pass

    def run_prep(self, map_dir, target, samples, abund_dir, **kwargs):
        bams = [join(map_dir, f'{target}_{s}.bam') for s in samples]
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

    def run_prep(self, assembly_dir, simulation_dir, abund_dir, samples, threads, **kwargs):
        # make sketches
        for s in samples:
            s = join(simulation_dir, s)
            cmd = f'{self.name} sketch 1 {s}_R1.fastq.gz -2 {s}_R2.fastq.gz -d {abund_dir}'
            run(cmd)

        # calculate coverage matrix
        contigs = join(assembly_dir, 'contigs.fasta')
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

    def run_main(self, simulation_dir, assembly_dir, abund_dir, target, map_sample, threads, **kwargs):
        # FORMAT FOR DIFFERENT BINNERS
        pass 

    def run_prep(self, simulation_dir, assembly_dir, abund_dir, map_sample, threads, **kwargs):
        contigs = join(assembly_dir, 'contigs.fasta')
        r1 = join(simulation_dir, f'{map_sample}_R1.fastq.gz')
        r2 = join(simulation_dir, f'{map_sample}_R2.fastq.gz')
        cmd = f'{self.execs} -t {threads} --aemb {contigs} {r1} {r2} > {abund_dir}/{map_sample}_coverage.tsv'
        run(cmd)

    def main_done(self, map_dir, target_sample, map_sample, **kwargs):
        # CHECK ALL FORMATS AVAILABLE
        pass

    def prep_done(self, abund_dir, map_sample, **kwargs):
        return not_empty(join(abund_dir, f'{map_sample}_coverage.tsv'))
