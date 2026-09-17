from abc import ABC, abstractmethod
import pandas as pd
from os import makedirs, rename
from os.path import join, exists, abspath, splitext, basename, dirname
from ..utils import not_empty, soft_link, rm_dir, run
from . import coverage as cv
import glob

def clear_prep(task_out_dir):
    rm_dir(join(task_out_dir, 'input'))

def clear_bin(task_out_dir):
    rm_dir(join(task_out_dir, 'output/bins'))

def task_done(bntsk, stage='all', report=False, only_not_done=False, clear=False):
    # TODO ADD QC STAGE, and BIN TABLE STAGE
    if stage == 'all' or stage == 'prep':
        task_out_dir = bntsk.task_out_dir
        b = getattr(bntsk.binner+'Binner')()
        prep_done =  b.prep_done(**bntsk.to_dict())
        prep_done &= exists(join(task_out_dir, 'prep_DONE'))
        if not prep_done and report:
            print('bin prep not DONE', bntsk.task_name, flush=True)
        if not prep_done and clear:
            clear_prep(task_out_dir)
            clear_bin(task_out_dir)
        if stage == 'prep':
            return prep_done

    if stage == 'all' or stage == 'bin':
        task_out_dir = bntsk.task_out_dir
        b = getattr(bn, (bntsk.binner + 'Binner') if not bntsk.is_refiner else (bntsk.refiner+'Refiner'))()
        bin_done =  b.bin_done(**bntsk.to_dict())
        bin_done &= exists(join(task_out_dir, 'bin_DONE'))
        if not bin_done and report:
            print('bin not DONE', bntsk.task_name, flush=True)
        if not bin_done and clear:
            clear_bin(task_out_dir)
        if stage == 'bin':
            return bin_done

    if stage == 'all' or stage == 'qc':
        task_out_dir = bntsk.task_out_dir
        b = getattr(qc, (bntsk.binner + 'Binner') if not bntsk.is_refiner else (bntsk.refiner+'Refiner'))()
        bin_done =  b.bin_done(**bntsk.to_dict())
        bin_done &= exists(join(task_out_dir, 'bin_DONE'))
        if not bin_done and report:
            print('bin not DONE', bntsk.task_name, flush=True)
        if not bin_done and clear:
            clear_bin(task_out_dir)
        if stage == 'bin':
            return bin_done

    return prep_done and bin_done

def bin_done(bin_dir):
    bins = glob.glob(f'{bin_dir}/*.fasta')
    return all(not_empty(e) for e in bins) and len(bins) > 0

def softlink_assembly_to_taskdir(asm_dir, sample, task_out_dir, prefix=None):
    contigs = abspath(join(asm_dir, f'contigs.fasta'))
    asm = join(task_out_dir, 'input', 'asm.fasta') if prefix is None else join(task_out_dir, 'input', f'{prefix}.fasta')
    soft_link(contigs, asm)
    return asm

def compute_idxstats(bam, idxstats):
    run(f'samtools idxstats {bam} > {idxstats}')

class Binner(ABC):
    name: str
    exec: str

    @abstractmethod
    def run_prep(self, asm_dir, samples, task_out_dir, **kwargs):
        ...

    @abstractmethod
    def run_main(self, task_out_dir, options, **kwargs):
        ...

    @abstractmethod
    def prep_done(self, task_out_dir, **kwargs) -> bool:
        ...

    @abstractmethod
    def main_done(self, task_out_dir, **kwargs) -> bool:
        ...

    @abstractmethod
    def bins_as_fasta(self, bin_dir) -> list:
        ...

class MaxBin2Binner(Binner):
    name='MaxBin2'
    exec='run_MaxBin.pl'
    #def run_main(self, task_out_dir, options, threads, **kwargs):
    #    makedirs(join(task_out_dir, f'output/bins'), exist_ok=True)
    #    contigs = join(task_out_dir, 'input', 'asm.fasta')
    #    assert(exists(contigs))
    #    counts = glob.glob(join(task_out_dir, 'input', '*.counts'))
    #    counts = ' '.join(f'-abund{i} {fl}' for i, fl in enumerate(counts,start=1))
    #    counts = counts.replace('-abund1', '-abund')
    #    cmd = f'{self.name} -thread {threads} {options} -contig {contigs} {counts} -out {task_out_dir}/output/bins/bin 2> {task_out_dir}/output/MaxBin2err.log'
    #    run(cmd)
    #    self.bins_as_fasta(f'{task_out_dir}/output/bins')

    def run_main(self, task_out_dir, options, threads, target, samples, coverage, cov_dir, **kwargs):
        apc = getattr(cv, f'{coverage}Coverage').abundance_precomputed()
        makedirs(join(task_out_dir, f'output/bins'), exist_ok=True)
        contigs = join(task_out_dir, 'input', 'asm.fasta')
        cov_dir = cov_dir if apc else join(task_out_dir, 'input')
        covs = [join(cov_dir, f'{s}_coverage.tsv') for s in samples]
        covs = ' '.join(f'-abund{i} {fl}' for i, fl in enumerate(covs,start=1))
        covs = covs.replace('-abund1', '-abund')
        cmd = f'{self.exec} -thread {threads} {options} -contig {contigs} {covs} -out {task_out_dir}/output/bins/bin 2> {task_out_dir}/output/MaxBin2err.log'
        run(cmd)
        self.bins_as_fasta(f'{task_out_dir}/output/bins')

    #def run_prep(self, asm_dir, map_dir, samples, task_out_dir, **kwargs):
    #    target_sample = samples[0]
    #    makedirs(join(task_out_dir, 'input'), exist_ok=True)
    #    # first get the contigs that will be binned (asm.fasta)
    #    softlink_assembly_to_taskdir(asm_dir, target_sample, task_out_dir)
    #    for s in samples:
    #        bam = join(map_dir, f'{target_sample}_{s}.bam')
    #        idxstats = join(task_out_dir, 'input', f'{s}.idxstats')
    #        # compute idx stats
    #        compute_idxstats(bam, idxstats)
    #        counts = pd.read_csv(idxstats, sep='\t', header=None, index_col=0)[2]
    #        counts = counts[:-1]  # Drop a trailing '*' from the idxstats file
    #        counts.to_csv(idxstats.replace('.idxstats','.counts'), sep='\t', header=None)
    
    def run_prep(self, asm_dir, cov_dir, samples, task_out_dir, coverage, **kwargs):
        apc = getattr(cv, f'{coverage}Coverage').abundance_precomputed()
        if apc: # nothing to do
            return
        target_sample = samples[0]
        makedirs(join(task_out_dir, 'input'), exist_ok=True)
        softlink_assembly_to_taskdir(asm_dir, target_sample, task_out_dir)
        bams = ' '.join([join(cov_dir, f'{e}.bam') for e in samples])
        cov = join(task_out_dir, 'input','coverage_mat_jgi.tsv')
        cv.jgi_summarize(bams, cov)
        cv.jgi_to_maxbin2(dirname(cov),cov)

    def bins_as_fasta(self, bin_dir):
        bins = glob.glob(join(bin_dir, '*.fasta'))
        return bins

    def prep_done(self, task_out_dir, samples, coverage, **kwargs):
        apc = getattr(cv, f'{coverage}Coverage').abundance_precomputed()
        if apc: return True
        fls = [join(task_out_dir, 'input', f'{s}_coverage.tsv') for s in samples]
        return all(not_empty(fl) for fl in fls)

    def main_done(self, task_out_dir, **kwargs):
        bin_dir = join(task_out_dir, 'output/bins')
        return bin_done(bin_dir)

    
class METABAT2Binner(Binner):
    name='METABAT2'
    exec='metabat2'

    def run_main(self,task_out_dir, options, threads, cov_dir, coverage, **kwargs):
        apc = getattr(cv, f'{coverage}Coverage').abundance_precomputed()
        cov_dir = cov_dir if apc else join(task_out_dir, 'input')
        makedirs(join(task_out_dir, f'output/bins'), exist_ok=True)
        contigs = join(task_out_dir, 'input', 'asm.fasta')
        cov = join(cov_dir, 'coverage_mat_jgi.tsv')
        cmd = f'{self.exec} --numThreads {threads} {options} --inFile {contigs}  --abdFile {cov} --outFile {task_out_dir}/output/bins/bin'
        run(cmd)
        self.bins_as_fasta(f'{task_out_dir}/output/bins')

    def run_prep(self, asm_dir, samples, task_out_dir, coverage, **kwargs):
        apc = getattr(cv, f'{coverage}Coverage').abundance_precomputed()
        if apc: # nothing to do
            return
        target_sample = samples[0]
        makedirs(join(task_out_dir, 'input'), exist_ok=True)
        softlink_assembly_to_taskdir(asm_dir, target_sample, task_out_dir)
        bams = ' '.join([join(asm_dir, f'{e}.bam') for e in samples])
        cv.jgi_summarize(bams, join(task_out_dir, 'input', 'coverage_mat_jgi.tsv'))

    def bins_as_fasta(self, bin_dir):
        bins = glob.glob(join(bin_dir, '*.fa'))
        if len(bins) != 0:
            for b in bins:
                rename(b, b.replace('.fa', '.fasta'))
        bins = glob.glob(join(bin_dir, '*.fasta'))
        return bins

    def prep_done(self, task_out_dir, coverage, **kwargs):
        apc = getattr(cv, f'{coverage}Coverage').abundance_precomputed()
        if apc: return True
        return not_empty(join(task_out_dir, 'input', 'coverage_mat_jgi.tsv'))

    def main_done(self, task_out_dir, **kwargs):
        bin_dir = join(task_out_dir, 'output/bins')
        return bin_done(bin_dir)


class VAMBBinner(Binner):
    name='VAMB'
    exec='/insomnia001/depts/pmg/users/ic2465/miniforge3/envs/vamb/bin/vamb'

    def run_main(self, task_out_dir, options, threads, cov_dir, coverage, **kwargs):
        contigs = join(task_out_dir, 'input', 'asm.fasta')
        makedirs(join(task_out_dir, 'output'),exist_ok=True)
        bin_dir = join(task_out_dir, 'output/bins')
        apc = getattr(cv, f'{coverage}Coverage').abundance_precomputed()
        abund = f'--abundance_tsv {join(cov_dir, 'coverage_mat.tsv')}' if apc else f'--bamdir {join(task_out_dir, 'input')}'
        cmd = f'{self.exec} bin default -p {threads} {options} --outdir {bin_dir} --fasta {contigs} {abund} --minfasta 100000 -o "" 2> {task_out_dir}/VAMB.errlog',
        run(cmd)
        self.bins_as_fasta(bin_dir)

    def run_prep(self, asm_dir, cov_dir, samples, task_out_dir, coverage, **kwargs):
        apc = getattr(cv, f'{coverage}Coverage').abundance_precomputed()
        if apc: # nothing to do
            return
        target_sample = samples[0]
        makedirs(join(task_out_dir, 'input'), exist_ok=True)
        softlink_assembly_to_taskdir(asm_dir, target_sample, task_out_dir)
        for s in samples:
            bam = join(cov_dir, f'{s}.bam')
            bamln = join(task_out_dir, 'input', f'{s}.bam')
            run(f'ln -s {bam} {bamln}')

    def bins_as_fasta(self, bin_dir):
        bins = glob.glob(join(bin_dir, 'bins/*.fna'))
        if len(bins) != 0:
            for b in bins:
                base = basename(b).replace('.fna','.fasta')
                bnew = join(bin_dir, f'bin.{base}')
                rename(b, bnew)
        bins = glob.glob(join(bin_dir, '*.fasta'))
        return bins

    def prep_done(self, task_out_dir, samples, coverage, **kwargs):
        apc = getattr(cv, f'{coverage}Coverage').abundance_precomputed()
        if apc: return True
        input_dir = join(task_out_dir, 'input')
        return not_empty(f'{input_dir}/asm.fasta') and all(not_empty(f'{input_dir}/{s}.bam') for s in samples)

    def main_done(self, task_out_dir, **kwargs):
        bin_dir = join(task_out_dir, 'output/bins')
        return bin_done(bin_dir)
    
class SemiBin2Binner(Binner):
    name='SemiBin2'
    exec = 'SemiBin2'

    def run_main(self, task_out_dir, options, target, samples, threads, cov_dir, coverage, **kwargs):
        apc = getattr(cv, f'{coverage}Coverage').abundance_precomputed()
        bin_dir = join(task_out_dir, 'output/bins')
        # single and multi-sample binning differ
        contigs = join(task_out_dir, 'input', 'asm.fasta')
        #if len(samples) == 1:
        #    bam = join(task_out_dir, 'input', f'{target_sample}_{target_sample}.bam')
        #    run(
        #        f'{self.exec} single_easy_bin --threads {threads} {options} --environment human_gut -i {contigs} -b {bam} ' + 
        #        f'-o {bin_dir} --compression none 2> {task_out_dir}/SemiBin2.errlog'
        #    )
        #else:
        if apc: 
            abund = '-a ' + ' '.join(f'{cov_dir}/{s}_coverage.tsv' for s in samples)
        else:
            abund = '-b ' + ' '.join(f'{task_out_dir}/input/{s}.bam' for s in samples)
        run(
            f'{self.exec} single_easy_bin --threads {threads} {options} -i {contigs} {abund} ' + 
            f'-o {bin_dir} --compression none 2> {task_out_dir}/SemiBin2.errlog'
        )
        self.bins_as_fasta(bin_dir)

    def run_prep(self, asm_dir, cov_dir, samples, target, task_out_dir, coverage, **kwargs):
        apc = getattr(cv, f'{coverage}coverage').abundance_precomputed()
        if apc: # nothing to do
            return
        makedirs(join(task_out_dir, 'input'), exist_ok=True)
        # single and multi-sample prep differ
        softlink_assembly_to_taskdir(asm_dir, target, task_out_dir)
        for sample in samples:
            bam = join(cov_dir, f'{sample}.bam')
            bamln = join(task_out_dir, 'input', f'{sample}.bam')
            run(f'ln -s {bam} {bamln}')

    def bins_as_fasta(self, bin_dir):
        bins = glob.glob(join(bin_dir, 'output_bins/*.fa'))
        if len(bins) != 0:
            for b in bins:
                base = basename(b).replace('.fa','.fasta').replace('SemiBin_', 'bin.')
                bnew = join(bin_dir, f'{base}')
                rename(b, bnew)
        bins = glob.glob(join(bin_dir, '*.fasta'))
        return bins

    def prep_done(self, task_out_dir, samples, coverage, **kwargs):
        apc = getattr(cv, f'{coverage}Coverage').abundance_precomputed()
        if apc: return True
        input_dir = join(task_out_dir, 'input')
        return not_empty(f'{input_dir}/asm.fasta') and all(not_empty(f'{input_dir}/{s}.bam') for s in samples)

    def main_done(self, task_out_dir, **kwargs):
        bin_dir = join(task_out_dir, 'output/bins')
        return bin_done(bin_dir)
    

class CONCOCTBinner(Binner):
    # Fits directly into current framework include
    name='CONCOCT'
    exec='/insomnia001/depts/pmg/users/ic2465/miniforge3/envs/concoct/bin/'

    def run_main(self, task_out_dir, options, threads, **kwargs):
        # run concoct main clustering
        input = join(task_out_dir, 'input')
        bin_dir = join(task_out_dir, 'output/bins')
        makedirs(bin_dir, exist_ok=True)
        run(
            f'{self.exec}/concoct --threads {threads} {options} --composition_file {input}/contig_10K.fa --coverage_file {input}/coverage_table.tsv -b {bin_dir} 2> {task_out_dir}/output/cncterr.log'
        )
        # merge subcontig clustering into orginal contig clustering
        run(
            f'{self.exec}/merge_cutup_clustering.py {bin_dir}/clustering_gt1000.csv > {bin_dir}/clustering_merged.csv 2> {bin_dir}/mergelog.txt'
        )
        # extract bins
        run(
            f'{self.exec}/extract_fasta_bins.py {input}/asm.fasta {bin_dir}/clustering_merged.csv --output_path {bin_dir}'
        )
        self.bins_as_fasta(bin_dir)

    def run_prep(self, asm_dir, cov_dir, target, samples, task_out_dir, **kwargs):
        makedirs(join(task_out_dir, 'input'), exist_ok=True)
        # first get the contigs that will be binned
        asm = softlink_assembly_to_taskdir(asm_dir, target, task_out_dir)
        # cut up the fasta into chunks
        run(
            f'{self.exec}/cut_up_fasta.py {asm} -c 10000 -o 0 --merge_last -b {task_out_dir}/input/contig_10K.bed ' +
            f'> {task_out_dir}/input/contig_10K.fa'
        )
        # make the coverage table
        bams = ' '.join([join(cov_dir, f'{target}_{s}.bam') for s in samples])
        run(
            f'{self.exec}/concoct_coverage_table.py {task_out_dir}/input/contig_10K.bed {bams} ' +
            f'> {task_out_dir}/input/coverage_table.tsv'
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

    def main_done(self, task_out_dir, **kwargs):
        bin_dir = join(task_out_dir, 'output/bins')
        return bin_done(bin_dir)


class COMEBinBinner:
    name='COMEBin'
    exec='comebin'
    size=1000
    def run_main(self, task_out_dir, options, threads, **kwargs):
        output = join(task_out_dir, 'output')
        input = join(task_out_dir, 'input')
        rm_dir(join(output, 'bins', 'comebin_res'))
        run(
            f'mamba run -n {self.exec} python -c "import torch; print(torch.cuda.is_available(),flush=True)"'
        )
        run(
            f'mamba run -n {self.exec} run_comebin.sh -t {threads} {options} -a {input}/asm_{self.size}.fa ' + 
            f'-p {input} -o {output}/bins -s main 2> {output}/COMEBin.errlog'
        )
        self.bins_as_fasta(f'{output}/bins')

    def run_prep(self, asm_dir, cov_dir, samples, task_out_dir, options, threads, **kwargs):
        input_dir = join(task_out_dir, 'input')
        output_dir = join(task_out_dir, 'output')
        makedirs(input_dir,exist_ok=True)
        makedirs(output_dir,exist_ok=True)
        makedirs(f'{output_dir}/bins', exist_ok=True)
        target_sample = samples[0]
        asm =  softlink_assembly_to_taskdir(asm_dir, target_sample, task_out_dir)
        # reduce to default (1000bp)
        run(f'mamba run -n {self.exec} /insomnia001/depts/pmg/users/ic2465/miniforge3/envs/comebin/bin/COMEBin/scripts/Filter_tooshort.py {asm} {self.size}')

        # make bed file for contigs > self.size
        df = pd.read_csv(f'{input_dir}/contig_length_filter{self.size}.txt', header=None, sep='\t')
        df.columns = ['contig_name', 'end']
        df['start'] = 0
        df[['contig_name', 'start', 'end']].to_csv(f'{input_dir}/contig_names{self.size}.bed',index=None, header=None, sep='\t')

        # extract bam for contigs > self.size
        for s in samples:
            run(
                f'samtools view -@ {threads} -b -L {input_dir}/contig_names{self.size}.bed {cov_dir}/{target_sample}_{s}.bam ' + 
                f'> {input_dir}/{target_sample}_{s}.{self.size}.bam'
            )

        # run data augmentation step
        run(
            f'mamba run -n {self.exec} run_comebin.sh -t {threads} {options} -a {input_dir}/asm_{self.size}.fa ' + 
            f'-p {input_dir} -o {input_dir}/tmp -s prep 2> {input_dir}/COMEBin.errlog'
        )
        bams = glob.glob(f'{input_dir}/*.bam')
        for b in bams:
            rm_file(b)
        move(f'{input}/tmp/data_augmentation', f'{output_dir}/bins/data_augmentation')

    def bins_as_fasta(self, bin_dir):
        bins = glob.glob(join(bin_dir, 'comebin_res/comebin_res_bins/*.fa'))
        if len(bins) != 0:
            for b in bins:
                base = basename(b).replace('.fa','.fasta')
                bnew = join(bin_dir, f'bin.{base}')
                rename(b, bnew)
        bins = glob.glob(join(bin_dir, '*.fasta'))
        return bins

    def prep_done(self, task_out_dir, **kwargs):
        input_dir = join(task_out_dir, 'output')
        means = glob.glob(f'{input_dir}/bins/data_augmentation/*_mean.tsv')
        vars = glob.glob(f'{input_dir}/bins/data_augmentation/*_var.tsv')
        return len(means) == len(vars) and all(not_empty(e) for e in means+vars) and len(means) > 0

    def main_done(self, task_out_dir, **kwargs):
        bin_dir = join(task_out_dir, 'output/bins')
        return bin_done(bin_dir)

class Refiner(ABC):
    name: str
    exec: str

    @abstractmethod
    def run_main(self, task_out_dir, options, **kwargs):
        ...

    @abstractmethod
    def main_done(self, task_out_dir, **kwargs) -> bool:
        ...

    @abstractmethod
    def bins_as_fasta(self, bin_dir) -> list:
        ...

class DAS_ToolRefiner(Refiner):
    name = 'DAS_Tool'
    exec = 'DAS_Tool'
    def run_main(self, task_out_dir, options, pipelines, ppln_bin_out,asm_dir, target_sample, threads=8, **kwargs):
        summaries = list()
        bin_dir = join(task_out_dir, f'output/bins')
        makedirs(bin_dir, exist_ok=True)
        input_dir = join(task_out_dir, 'input') 
        makedirs(input_dir, exist_ok=True)

        for i, bin_out in enumerate(ppln_bin_out):
            bins = glob.glob(join(bin_out, 'output/bins/*.fasta'))
            das_summary = join(input_dir, f'ppln_{i}_das_summary.tsv')
            hdrs=list()
            for bin in bins:
                with open(bin) as fin:
                    bin_name = splitext(basename(bin))[0]
                    hdrs += [(e.strip()[1:], bin_name) for e in fin if e.startswith('>')]
            pd.DataFrame(hdrs).to_csv(das_summary, header=None, index=None, sep='\t')
            summaries.append(das_summary)
        contigs = softlink_assembly_to_taskdir(asm_dir, target_sample, task_out_dir)
        contig2bin = ','.join(summaries)
        run(
            f'DAS_Tool -t {threads} {options} -i {contig2bin} -c {contigs} -o {bin_dir} --write_bin_evals --write_bins --write_unbinned'
        )
        self.bins_as_fasta(bin_dir)
            
    def bins_as_fasta(self, bin_dir):
        bins = glob.glob(bin_dir+ '_DASTool_bins/*.fa')
        bins = [e for e in bins if 'unbinned.fa' not in e]
        records = list()
        if len(bins) != 0:
            for i, b in enumerate(bins):
                base = basename(b)
                bnew = join(bin_dir, f'bin.{i}.fasta')
                records.append((base, f'bin.{i}.fasta'))
                rename(b, bnew)
            pd.DataFrame(records, columns=['das_name', 'new_name']).to_csv(join(bin_dir, f'bin_renaming_mapping.csv'),index=None)
        bins = glob.glob(join(bin_dir, 'bin*.fasta'))
        return bins

    @staticmethod
    def main_done(task_out_dir, **kwargs):
        bin_dir = join(task_out_dir, 'output/bins')
        return bin_done(bin_dir)
