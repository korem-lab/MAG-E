from abc import ABC, abstractmethod
import pandas as pd
from os import makedirs, rename
from os.path import join, exists, abspath, splitext, basename, dirname
from subprocess import run
from ..utils import not_empty
import glob

def jgi_summarize(bams, filename):
    log = splitext(filename)[0] + '.log'
    if exists(filename):
        return
    run(f'jgi_summarize_bam_contig_depths --outputDepth {filename} {bams} 2> {log}',shell=True)

def bin_done(bin_dir):
    bins = glob.glob(f'{bin_dir}/*.fasta')
    return len(bins) > 0 and all(not_empty(e) for e in bins) 

def softlink_assembly_to_taskdir(asm_dir, sample, task_out_dir, prefix=None):
    contigs = abspath(join(asm_dir, f'{sample}.contigs.fasta'))
    asm = join(task_out_dir, 'input', 'asm.fasta') if prefix is None else join(task_out_dir, 'input', f'{prefix}.fasta')
    run(f'ln -f -s {contigs} {asm}', shell=True)
    return asm

def compute_idxstats(bam, idxstats):
    run(f'samtools idxstats {bam} > {idxstats}',shell=True)

class Binner(ABC):
    name: str
    execs: str

    @abstractmethod
    def run_prep(self, asm_dir, bam_dir, samples, task_out_dir, **kwargs):
        ...

    @abstractmethod
    def run_binning(self, task_out_dir, options, **kwargs):
        ...

    @abstractmethod
    def construct_binning_table(self, task_name, sample, task_out_dir, ground_truth):
        ...

    @abstractmethod
    def prep_done(self, task_out_dir, **kwargs) -> bool:
        ...

    @abstractmethod
    def bin_done(self, task_out_dir, **kwargs) -> bool:
        ...

    @abstractmethod
    def bin_as_fasta(self, bin_dir) -> list:
        ...

class MaxBin2Binner(Binner):
    name='MaxBin2'
    exec='run_MaxBin.pl'
    def run_binning(self, task_out_dir, binning_options, threads, **kwargs):
        makedirs(join(task_out_dir, f'output/bins'), exist_ok=True)
        contigs = join(task_out_dir, 'input', 'asm.fasta')
        assert(exists(contigs))
        counts = glob.glob(join(task_out_dir, 'input', '*.counts'))
        counts = ' '.join(f'-abund{i} {fl}' for i, fl in enumerate(counts,start=1))
        counts = counts.replace('-abund1', '-abund')
        cmd = f'{self.name} -thread {threads} {binning_options} -contig {contigs} {counts} -out {task_out_dir}/output/bins/bin 2> {task_out_dir}/output/MaxBin2err.log'
        run(cmd, shell=True)
        self.bins_as_fasta(f'{task_out_dir}/output/bins')

    def run_prep(self, asm_dir, samples, task_out_dir, **kwargs):
        target_sample = samples[0]
        makedirs(join(task_out_dir, 'input'), exist_ok=True)
        # first get the contigs that will be binned (asm.fasta)
        softlink_assembly_to_taskdir(asm_dir, target_sample, task_out_dir)
        for s in samples:
            bam = join(asm_dir, f'{target_sample}_{s}.bam')
            idxstats = join(task_out_dir, 'input', f'{s}.idxstats')
            # compute idx stats
            compute_idxstats(bam, idxstats)
            counts = pd.read_csv(idxstats, sep='\t', header=None, index_col=0)[2]
            counts = counts[:-1]  # Drop a trailing '*' from the idxstats file
            counts.to_csv(idxstats.replace('.idxstats','.counts'), sep='\t', header=None)
    
    def bins_as_fasta(self, bin_dir):
        bins = glob.glob(join(bin_dir, '*.fasta'))
        return bins

    def prep_done(self, task_out_dir, samples, **kwargs):
        return all(not_empty(join(task_out_dir, 'input', f'{s}.counts')) for s in samples)

    def bin_done(self, task_out_dir, **kwargs):
        bin_dir = join(task_out_dir, 'output/bins')
        return bin_done(bin_dir)

    
class METABAT2Binner(Binner):
    name='METABAT2'
    exec='metabat2'

    def run_binning(self,task_out_dir, binning_options, threads, **kwargs):
        makedirs(join(task_out_dir, f'output/bins'), exist_ok=True)
        contigs = join(task_out_dir, 'input', 'asm.fasta')
        count_mat = join(task_out_dir, 'input', 'count_mat.tsv')
        cmd = f'{self.exec} --numThreads {threads} {binning_options} --inFile {contigs}  --abdFile {count_mat} --outFile {task_out_dir}/output/{METABAT2Binner.name}_bins/bin'
        run(cmd, shell=True)
        self.bins_as_fasta(f'{task_out_dir}/output/bins')

    def run_prep(self, asm_dir, bam_dir, samples, task_out_dir, **kwargs):
        target_sample = samples[0]
        makedirs(join(task_out_dir, 'input'), exist_ok=True)
        softlink_assembly_to_taskdir(asm_dir, target_sample, task_out_dir)
        bams = ' '.join([join(bam_dir, f'{target_sample}_{e}.bam') for e in samples])
        jgi_summarize(bams, join(task_out_dir, 'input', 'count_mat.tsv'))

    def bins_as_fasta(self, bin_dir):
        bins = glob.glob(join(bin_dir, '*.fa'))
        if len(bins) != 0:
            for b in bins:
                rename(b, b.replace('.fa', '.fasta'))
        bins = glob.glob(join(bin_dir, '*.fasta'))
        return bins

    def prep_done(self, task_out_dir, **kwargs):
        return not_empty(join(task_out_dir, 'input', 'count_mat.tsv'))

    def bin_done(self, task_out_dir, **kwargs):
        bin_dir = join(task_out_dir, 'output/bins')
        return bin_done(bin_dir)


class VAMBBinner(Binner):
    name='VAMB'
    execs='vamb'

    def run_binning(self, task_out_dir, options, **kwargs):
        input = join(task_out_dir, 'input')
        bin_dir = join(task_out_dir, 'output/bins')
        makedirs(join(task_out_dir, 'output'),exist_ok=True)
        cmd = f'{self.execs} bin default {options} --outdir {bin_dir} --fasta {input}/asm.fasta --bamdir {input} --minfasta 100000 -o "" 2> {task_out_dir}/VAMB.errlog',
        run(cmd, shell=True)
        self.bins_as_fasta(bin_dir)

    def run_prep(self, asm_dir, bam_dir, samples, task_out_dir, **kwargs):
        target_sample = samples[0]
        makedirs(join(task_out_dir, 'input'), exist_ok=True)
        softlink_assembly_to_taskdir(asm_dir, target_sample, task_out_dir)
        for s in samples:
            bam = join(bam_dir, f'{target_sample}_{s}.bam')
            bamln = join(task_out_dir, 'input', f'{target_sample}_{s}.bam')
            run(f'ln -s {bam} {bamln}',shell=True)

    def bins_as_fasta(self, bin_dir):
        bins = glob.glob(join(bin_dir, 'bins/*.fna'))
        if len(bins) != 0:
            for b in bins:
                base = basename(b).replace('.fna','.fasta')
                bnew = join(bin_dir, f'bin.{base}')
                rename(b, bnew)
        bins = glob.glob(join(bin_dir, '*.fasta'))
        return bins

    def prep_done(self, task_out_dir, samples, **kwargs):
        input_dir = join(task_out_dir, 'input')
        ts = samples[0]
        return not_empty(f'{input_dir}/asm.fasta') and all(not_empty(f'{input_dir}/{ts}_{s}.bam') for s in samples)

    def bin_done(self, task_out_dir, **kwargs):
        bin_dir = join(task_out_dir, 'output/bins')
        return bin_done(bin_dir)
    
class SemiBin2Binner(Binner):
    name='SemiBin2'
    exec = 'SemiBin2'

    def run_binning(self, task_out_dir, options, samples, **kwargs):
        bin_dir = join(task_out_dir, 'output/bins')
        # single and multi-sample binning differ
        contigs = join(task_out_dir, 'input', 'asm.fasta')
        target_sample = samples[0]
        if len(samples) == 1:
            bam = join(task_out_dir, 'input', f'{target_sample}_{target_sample}.bam')
            run(
                f'{self.exec} single_easy_bin {options} --environment human_gut -i {contigs} -b {bam} ' + 
                f'-o {bin_dir} --compression none 2> {task_out_dir}/SemiBin2.errlog',
                shell=True
            )
        else:
            bams = ' '.join(f'{task_out_dir}/input/{target_sample}_{s}.bam' for s in samples)
            run(
                f'{self.exec} single_easy_bin {options} -i {contigs} -b {bams} ' + 
                f'-o {bin_dir} --compression none 2> {task_out_dir}/SemiBin2.errlog',
                shell=True
            )
        SemiBin2Binner.bins_as_fasta(bin_dir)

    def run_prep(self, asm_dir, bam_dir, samples, task_out_dir, **kwargs):
        makedirs(join(task_out_dir, 'input'), exist_ok=True)
        # single and multi-sample prep differ
        if len(samples) == 1:
            target_sample = samples[0]
            bam = join(bam_dir, f'{target_sample}_{target_sample}.bam')
            bamln = join(task_out_dir, 'input', f'{target_sample}_{target_sample}.bam')
            run(f'ln -s {bam} {bamln}',shell=True)
            # first get the contigs that will be binned
            softlink_assembly_to_taskdir(asm_dir, target_sample, task_out_dir)
        else:
            target_sample = samples[0]
            softlink_assembly_to_taskdir(asm_dir, target_sample, task_out_dir)
            for sample in samples:
                bam = join(bam_dir, f'{target_sample}_{sample}.bam')
                bamln = join(task_out_dir, 'input', f'{target_sample}_{sample}.bam')
                run(f'ln -s {bam} {bamln}',shell=True)


    def bins_as_fasta(self, bin_dir):
        bins = glob.glob(join(bin_dir, 'output_bins/*.fa'))
        if len(bins) != 0:
            for b in bins:
                base = basename(b).replace('.fa','.fasta').replace('SemiBin_', 'bin.')
                bnew = join(bin_dir, f'{base}')
                rename(b, bnew)
        bins = glob.glob(join(bin_dir, '*.fasta'))
        return bins

    def prep_done(self, task_out_dir, samples, **kwargs):
        input_dir = join(task_out_dir, 'input')
        return not_empty(f'{input_dir}/asm.fasta') and all(not_empty(f'{input_dir}/{samples[0]}_{s}.bam') for s in samples)

    def bin_done(self, task_out_dir, **kwargs):
        bin_dir = join(task_out_dir, 'output/bins')
        return bin_done(bin_dir)
    

class CONCOCTBinner(Binner):
    # Fits directly into current framework include
    name='CONCOCT'
    exec='/concoct/dir'

    def run_binning(self, task_out_dir, options, **kwargs):
        # run concoct main clustering
        input = join(task_out_dir, 'input')
        bin_dir = join(task_out_dir, 'output/bins')
        makedirs(bin_dir, exist_ok=True)
        run(
            f'{self.execs}/concoct {options} --composition_file {input}/contig_10K.fa --coverage_file {input}/coverage_table.tsv -b {bin_dir} 2> {task_out_dir}/output/cncterr.log',
            shell=True
        )
        # merge subcontig clustering into orginal contig clustering
        run(
            f'{self.execs}/merge_cutup_clustering.py {bin_dir}/clustering_gt1000.csv > {bin_dir}/clustering_merged.csv 2> {bin_dir}/mergelog.txt',
            shell=True
        )
        # extract bins
        run(
            f'{self.execs}/extract_fasta_bins.py {input}/asm.fasta {bin_dir}/clustering_merged.csv --output_path {bin_dir}',
            shell=True
        )
        self.bins_as_fasta(bin_dir)

    def run_prep(self, asm_dir, bam_dir, samples, task_out_dir, **kwargs):
        target_sample = samples[0]
        makedirs(join(task_out_dir, 'input'), exist_ok=True)
        # first get the contigs that will be binned
        asm = softlink_assembly_to_taskdir(asm_dir, target_sample, task_out_dir)
        # cut up the fasta into chunks
        run(
            f'{self.execs}/cut_up_fasta.py {asm} -c 10000 -o 0 --merge_last -b {task_out_dir}/input/contig_10K.bed ' +
            f'> {task_out_dir}/input/contig_10K.fa',
            shell=True
        )
        # make the coverage table
        bams = ' '.join([join(bam_dir, f'{target_sample}_{s}.bam') for s in samples])
        run(
            f'{self.execs}/concoct_coverage_table.py {task_out_dir}/input/contig_10K.bed {bams} ' +
            f'> {task_out_dir}/input/coverage_table.tsv',
            shell=True
        )

    def bins_as_fasta(self, bin_dir):
        bins = glob.glob(join(bin_dir, '*.fa'))
        if len(bins) != 0:
            for b in bins:
                base = basename(b).replace('.fa','.fasta')
                bnew = join(dirname(b), f'bin.{base}')
                rename(b, bnew)
        bins = glob.glob(join(bin_dir, '*.fasta'))
        return bins

    def prep_done(self, task_out_dir, **kwargs):
        input_dir = join(task_out_dir, 'input')
        return not_empty(f'{input_dir}/contig_10K.bed') and not_empty(f'{input_dir}/contig_10K.fa') and not_empty(f'{input_dir}/coverage_table.tsv')

    def bin_done(self, task_out_dir, **kwargs):
        bin_dir = join(task_out_dir, 'output/bins')
        return bin_done(bin_dir)


class COMEBinBinner:
    name='COMEBin'
    execs='comebin'
    size=1000
    def run_binning(self, task_out_dir, options, **kwargs):
        output = join(task_out_dir, 'output')
        input = join(task_out_dir, 'input')
        makedirs(output, exist_ok=True)

        run(
            f'mamba run -n {self.execs} run_comebin.sh {options} -a {input}/asm_{self.execs}.fa ' + 
            f'-p {input} -o {output}/bins 2> {output}/COMEBin.errlog', shell=True
        )
        self.bins_as_fasta(f'{output}/bins')

    def run_prep(self, asm_dir, bam_dir, samples, task_out_dir, **kwargs):
        input_dir = join(task_out_dir, 'input')
        makedirs(input_dir,exist_ok=True)
        target_sample = samples[0]
        asm =  softlink_assembly_to_taskdir(asm_dir, target_sample, task_out_dir)
        # reduce to default (1000bp)
        run(f'mamba run -n {self.execs} Filter_tooshort.py {asm} {self.size}',shell=True)

        # make bed file for contigs > self.size
        df = pd.read_csv(f'{input_dir}/contig_length_filter{self.size}.txt', header=None, sep='\t')
        df.columns = ['contig_name', 'end']
        df['start'] = 0
        df[['contig_name', 'start', 'end']].to_csv(f'{input_dir}/contig_names{self.size}.bed',index=None, header=None, sep='\t')

        # extract bam for contigs > self.size
        for s in samples:
            run(
                f'samtools view -@ 4 -b -L {input_dir}/contig_names{self.size}.bed {bam_dir}/{target_sample}_{s}.bam ' + 
                f'> {input_dir}/{target_sample}_{s}.{self.size}.bam', shell=True
            )
        
    def bins_as_fasta(self, bin_dir):
        bins = glob.glob(join(bin_dir, 'comebin_res/comebin_res_bins/*.fa'))
        if len(bins) != 0:
            for b in bins:
                base = basename(b).replace('.fa','.fasta')
                bnew = join(bin_dir, f'bin.{base}')
                rename(b, bnew)
        bins = glob.glob(join(bin_dir, '*.fasta'))
        return bins

    def prep_done(self, task_out_dir, samples, **kwargs):
        input_dir = join(task_out_dir, 'input')
        ts = samples[0]
        return not_empty(f'{input_dir}/asm_{self.size}.fa') and all(not_empty(f'{input_dir}/{ts}_{s}.{self.size}.bam') for s in samples)

    def bin_done(self, task_out_dir, **kwargs):
        bin_dir = join(task_out_dir, 'output/bins')
        return bin_done(bin_dir)

class DAS_ToolRefiner:
    pass