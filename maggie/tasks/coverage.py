from subprocess import run, DEVNULL
from glob import glob
import pandas as pd
from abc import ABC, abstractmethod
from os.path import splitext, join, dirname, basename
from os import makedirs
from ..utils import rm_dir, rm_file, not_empty, run, get_temp_local

def sort_bam(bam, coordinate=True, tmp_pref='tmp', threads=8):
    sort_coordinate = '' if coordinate else ' -n '
    dirname  = get_temp_local()
    makedirs(dirname, exist_ok=True)
    run(f'samtools sort -@ {threads} -m 2G -T {dirname} {sort_coordinate} {bam} > {bam}.{tmp_pref}')
    run(f'mv {bam}.{tmp_pref} {bam}')

def index_bam(bam):
    bai = bam.replace('.bam', '.bai')
    run(f'samtools index {bam} {bai}')

def jgi_summarize(bams, filename):
    log = splitext(filename)[0] + '.log'
    run(f'jgi_summarize_bam_contig_depths --outputDepth {filename} {bams} 2> {log}')

def jgi_to_maxbin2(dir, coverage, pref=''):
    df = pd.read_csv(coverage, delimiter='\t')
    for s in df.columns[3::2]:
        df[['contigName',s]].to_csv(
            join(dir, f'{pref}{s.replace('.bam', '')}_coverage.tsv'),
            index=None, header=None, sep='\t'
        )

def jgi_to_vamb(dir, coverage,pref=''):
    df = pd.read_csv(coverage, delimiter='\t')
    cov_cols = df.columns[3::2]
    cov_df = df[['contigName'] + list(cov_cols)]
    cov_df.columns = [e.replace('.bam', '') for e in cov_df.columns]
    cov_df.rename({'contigName':'contigname'},axis=1,inplace=True)
    cov_df.to_csv(join(dir, f'{pref}coverage_mat.tsv'), sep='\t', index=None)

def aembs_to_vamb(dir, covs, pref=''):
    cov_mat = pd.concat(
        [pd.read_csv(e, index_col=0, delimiter='\t', header=None) for e in covs], axis=1
    )
    cov_mat.columns = [basename(e).split('_')[0] for e in covs]
    cov_mat.index.name = 'contigname'
    cov_mat.to_csv(join(dir, f'{pref}coverage_mat.tsv'), sep='\t')

def make_formats(bams, jgi, pref):
    jgi_summarize(bams, jgi)
    jgi_to_vamb(dirname(jgi), jgi, pref)
    jgi_to_maxbin2(dirname(jgi), jgi, pref)

class Coverage(ABC):
    name: str
    exec: str

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

    def abundance_precomputed():
        return False

    def run_main(self, simulation_dir, cov_dir, cov_sample, threads, options, **kwargs):
        r1 = join(simulation_dir, f'{cov_sample}_R1.fastq.gz')
        r2 = join(simulation_dir, f'{cov_sample}_R2.fastq.gz')
        bam = join(cov_dir, f'{cov_sample}.bam')
        cmd = f'bowtie2 -p {threads} {options} -x {cov_dir}/idx -1 {r1} -2 {r2} | samtools view -bS - > {bam}'
        run(cmd)
        sort_bam(bam, threads)
        index_bam(bam)

    def run_prep(self, cov_dir, asm_dir, threads, **kwargs):
        rm_dir(cov_dir, remake=True)
        cmd = f'bowtie2-build --threads {threads} {asm_dir}/contigs.fasta {cov_dir}/idx'
        run(cmd)
    
    def main_done(self, cov_dir, cov_sample, **kwargs):
        bam = join(cov_dir, f'{cov_sample}.bam')
        return not_empty(bam)

    def prep_done(self, cov_dir, **kwargs):
        # check the assembly is present
        files = glob(f'{cov_dir}/*.bt2') + glob(f'{cov_dir}/*.bt2l')
        return all(not_empty(f) for f in files) and files

class bowtie2JGICoverage(Coverage):

    def abundance_precomputed():
        return True

    def run_main(self, cov_dir, target, **kwargs):
        # get to the coverage file 
        bams = glob(join(cov_dir, '*.bam'))
        bams = ' '.join(bams)
        jgi = join(cov_dir, 'coverage_mat_jgi.tsv')
        make_formats(bams, jgi, '')

        # single mode will be collected in same coverage_task directory
        bams = join(cov_dir, f'{target}.bam')
        jgi = join(cov_dir, 'target_coverage_mat_jgi.tsv')
        make_formats(bams, jgi, 'target_')

        # cleanup
        for f in glob(join(cov_dir, '*.bam')):
            rm_file(f)
        for f in glob(join(cov_dir, '*.bai')):
            rm_file(f)

    def run_prep(self, cov_dir, simulation_dir, samples, asm_dir, threads, options, **kwargs):
        # remove old files
        rm_dir(cov_dir, remake=True)

        # index
        cmd = f'bowtie2-build --threads {threads} {asm_dir}/contigs.fasta {cov_dir}/idx'
        run(cmd)

        # map
        for s in samples:
            r1 = join(simulation_dir, f'{s}_R1.fastq.gz')
            r2 = join(simulation_dir, f'{s}_R2.fastq.gz')
            bam = join(cov_dir, f'{s}.bam')
            cmd = f'bowtie2 -p {threads} -u 100000 {options} -x {cov_dir}/idx -1 {r1} -2 {r2} | samtools view -bS - > {bam}'
            run(cmd)
            sort_bam(bam, threads)
            index_bam(bam)

    def main_done(self, cov_dir, samples, target, **kwargs):
        fls = [join(cov_dir, 'coverage_mat_jgi.tsv'), join(cov_dir, 'coverage_mat.tsv')] + \
              [join(cov_dir, 'target_coverage_mat_jgi.tsv'), join(cov_dir, 'target_coverage_mat.tsv')] + \
            [join(cov_dir, f'{s}_coverage.tsv') for s in samples] + [join(cov_dir, f'target_{target}_coverage.tsv')]
        return all(not_empty(e) for e in fls)

    def prep_done(self, cov_dir, samples, **kwargs):
        return all(not_empty(join(cov_dir, f'{s}.bam')) for s in samples)

class bwamem2JGICoverage(Coverage):

    def abundance_precomputed():
        return True

    def run_main(self, cov_dir, target, **kwargs):
        # get to the coverage file 
        bams = glob(join(cov_dir, '*.bam'))
        bams = ' '.join(bams)
        jgi = join(cov_dir, 'coverage_mat_jgi.tsv')
        make_formats(bams, jgi,'')

        # single mode will be collected in same coverage_task directory
        bams = join(cov_dir, f'{target}.bam')
        jgi = join(cov_dir, 'target_coverage_mat_jgi.tsv')
        make_formats(bams, jgi, 'target_')

        # cleanup
        for f in glob(join(cov_dir, '*.bam')):
            rm_file(f)
        for f in glob(join(cov_dir, '*.bai')):
            rm_file(f)

    def run_prep(self, cov_dir, simulation_dir, samples, asm_dir, threads, options, **kwargs):
        # remove old files
        rm_dir(cov_dir, remake=True)

        # index
        cmd = f'bwa-mem2 index -p {cov_dir}/idx {asm_dir}/contigs.fasta '
        run(cmd)

        # map
        for s in samples:
            r1 = join(simulation_dir, f'{s}_R1.fastq.gz')
            r2 = join(simulation_dir, f'{s}_R2.fastq.gz')
            bam = join(cov_dir, f'{s}.bam')
            cmd = f'bwa-mem2 mem -t {threads} {options}  {cov_dir}/idx {r1} {r2} | samtools view -bS - > {bam}'
            run(cmd)
            sort_bam(bam, threads)
            index_bam(bam)

    def main_done(self, cov_dir, samples, target, **kwargs):
        fls = [join(cov_dir, 'coverage_mat_jgi.tsv'), join(cov_dir, 'coverage_mat.tsv')] + \
              [join(cov_dir, 'target_coverage_mat_jgi.tsv'), join(cov_dir, 'target_coverage_mat.tsv')] + \
            [join(cov_dir, f'{s}_coverage.tsv') for s in samples] + [join(cov_dir, f'target_{target}_coverage.tsv')]
        return all(not_empty(e) for e in fls)

    def prep_done(self, cov_dir, samples, **kwargs):
        return all(not_empty(join(cov_dir, f'{s}.bam')) for s in samples)

class fairyCoverage(Coverage):
    name= 'fairy'
    exec= 'fairy'

    def abundance_precomputed():
        return True

    def run_main(self, asm_dir, cov_dir, threads, target, **kwargs):
        # calculate coverage matrix
        contigs = join(asm_dir, 'contigs.fasta')
        cmd = f'fairy coverage {cov_dir}/*.bcsp {contigs} -t {threads} -o {cov_dir}/coverage_mat_jgi.tsv'
        run(cmd)
        jgi_to_maxbin2(cov_dir, join(cov_dir, 'coverage_mat_jgi.tsv'))
        jgi_to_vamb(cov_dir, join(cov_dir, 'coverage_mat_jgi.tsv'))

        # single mode
        contigs = join(asm_dir, 'contigs.fasta')
        cmd = f'fairy coverage {cov_dir}/{target}.paired.bcsp {contigs} -t {threads} -o {cov_dir}/target_coverage_mat_jgi.tsv'
        run(cmd)
        jgi_to_maxbin2(cov_dir, join(cov_dir, 'target_coverage_mat_jgi.tsv'), 'target_')
        jgi_to_vamb(cov_dir, join(cov_dir, 'target_coverage_mat_jgi.tsv'), 'target_')

    def run_prep(self, simulation_dir, cov_dir, samples, threads, **kwargs):
        # make sketches
        for s in samples:
            r = join(simulation_dir, s)
            cmd = f'fairy sketch -t {threads} -1 {r}_R1.fastq.gz -2 {r}_R2.fastq.gz -S {s} -d {cov_dir}'
            run(cmd)

    def main_done(self, cov_dir, samples, target, **kwargs):
        fls = [join(cov_dir, 'coverage_mat_jgi.tsv'), join(cov_dir, 'coverage_mat.tsv')] + \
              [join(cov_dir, 'target_coverage_mat_jgi.tsv'), join(cov_dir, 'target_coverage_mat.tsv')] + \
            [join(cov_dir, f'{s}_coverage.tsv') for s in samples] + [join(cov_dir, f'target_{target}_coverage.tsv')]
        return all(not_empty(e) for e in fls)

    def prep_done(self, cov_dir, samples, **kwargs):
        return all(not_empty(join(cov_dir, f'{s}.paired.bcsp')) for s in samples)

class aembCoverage(Coverage):
    name= 'strobealign'
    exec= 'strobealign'

    def abundance_precomputed():
        return True

    def run_main(self, **kwargs):
        pass 

    def run_prep(self, simulation_dir, asm_dir, cov_dir, samples, target, threads, **kwargs):
        contigs = join(asm_dir, 'contigs.fasta')
        for s in samples:
            r1 = join(simulation_dir, f'{s}_R1.fastq.gz')
            r2 = join(simulation_dir, f'{s}_R2.fastq.gz')
            cmd = f'{self.exec} -t {threads} --aemb {contigs} {r1} {r2} -o {cov_dir}/{s}_coverage.tsv'
            run(cmd)
        covs = [join(cov_dir, f'{s}_coverage.tsv') for s in samples]
        aembs_to_vamb(cov_dir, covs)

        # single mode
        r1 = join(simulation_dir, f'{target}_R1.fastq.gz')
        r2 = join(simulation_dir, f'{target}_R2.fastq.gz')
        cmd = f'{self.exec} -t {threads} --aemb {contigs} {r1} {r2} -o {cov_dir}/target_{target}_coverage.tsv'
        run(cmd)
        covs = [join(cov_dir, f'{s}_coverage.tsv') for s in [target]]
        aembs_to_vamb(cov_dir, covs,'target_')


    def main_done(self, **kwargs):
        return True

    def prep_done(self, cov_dir, samples, target, **kwargs):
        fls = [join(cov_dir, 'coverage_mat.tsv'),join(cov_dir, 'target_coverage_mat.tsv')] + \
            [join(cov_dir, f'{s}_coverage.tsv') for s in samples] + [join(cov_dir, f'target_{target}_coverage.tsv')]
        return all(not_empty(e) for e in fls)
