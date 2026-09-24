from abc import ABC, abstractmethod
import pandas as pd
from os import makedirs, rename
from os.path import join, exists, abspath, splitext, basename, dirname
from ..utils import not_empty, soft_link, rm_dir, run
from . import coverage as cv
import glob

def clear_prep(binner_dir):
    rm_dir(join(binner_dir, 'input'))

def clear_bin(binner_dir):
    rm_dir(join(binner_dir, 'output/bins'))

def task_done(bntsk, stage='all', report=False, only_not_done=False, clear=False):
    # TODO ADD QC STAGE, and BIN TABLE STAGE
    if stage == 'all' or stage == 'prep':
        binner_dir = bntsk.binner_dir
        b = getattr(bntsk.binner+'Binner')()
        prep_done =  b.prep_done(**bntsk.to_dict())
        prep_done &= exists(join(binner_dir, 'prep_DONE'))
        if not prep_done and report:
            print('bin prep not DONE', bntsk.task_hash, flush=True)
        if not prep_done and clear:
            clear_prep(binner_dir)
            clear_bin(binner_dir)
        if stage == 'prep':
            return prep_done

    if stage == 'all' or stage == 'bin':
        binner_dir = bntsk.binner_dir
        b = getattr(bn, (bntsk.binner + 'Binner') if not bntsk.is_refiner else (bntsk.refiner+'Refiner'))()
        bin_done =  b.bin_done(**bntsk.to_dict())
        bin_done &= exists(join(binner_dir, 'bin_DONE'))
        if not bin_done and report:
            print('bin not DONE', bntsk.task_hash, flush=True)
        if not bin_done and clear:
            clear_bin(binner_dir)
        if stage == 'bin':
            return bin_done

    if stage == 'all' or stage == 'qc':
        binner_dir = bntsk.binner_dir
        b = getattr(qc, (bntsk.binner + 'Binner') if not bntsk.is_refiner else (bntsk.refiner+'Refiner'))()
        bin_done =  b.bin_done(**bntsk.to_dict())
        bin_done &= exists(join(binner_dir, 'bin_DONE'))
        if not bin_done and report:
            print('bin not DONE', bntsk.task_hash, flush=True)
        if not bin_done and clear:
            clear_bin(binner_dir)
        if stage == 'bin':
            return bin_done

    return prep_done and bin_done

def bin_done(bin_dir):
    bins = glob.glob(f'{bin_dir}/*.fasta')
    return all(not_empty(e) for e in bins) and len(bins) > 0

def softlink_assembly_to_taskdir(assembly_dir, sample, binner_dir, prefix=None):
    contigs = abspath(join(assembly_dir, f'contigs.fasta'))
    asm = join(binner_dir, 'input', 'asm.fasta') if prefix is None else join(binner_dir, 'input', f'{prefix}.fasta')
    soft_link(contigs, asm)
    return asm

def compute_idxstats(bam, idxstats):
    run(f'samtools idxstats {bam} > {idxstats}')

class Binner(ABC):
    name: str
    exec: str

    @abstractmethod
    def run_prep(self, assembly_dir, samples, binner_dir, **kwargs):
        ...

    @abstractmethod
    def run_main(self, binner_dir, options, **kwargs):
        ...

    @abstractmethod
    def prep_done(self, binner_dir, **kwargs) -> bool:
        ...

    @abstractmethod
    def main_done(self, binner_dir, **kwargs) -> bool:
        ...

    @abstractmethod
    def bins_as_fasta(self, bin_dir) -> list:
        ...

class MaxBin2Binner(Binner):
    name='MaxBin2'
    exec='run_MaxBin.pl'
    #def run_main(self, binner_dir, options, threads, **kwargs):
    #    makedirs(join(binner_dir, f'output/bins'), exist_ok=True)
    #    contigs = join(binner_dir, 'input', 'asm.fasta')
    #    assert(exists(contigs))
    #    counts = glob.glob(join(binner_dir, 'input', '*.counts'))
    #    counts = ' '.join(f'-abund{i} {fl}' for i, fl in enumerate(counts,start=1))
    #    counts = counts.replace('-abund1', '-abund')
    #    cmd = f'{self.name} -thread {threads} {options} -contig {contigs} {counts} -out {binner_dir}/output/bins/bin 2> {binner_dir}/output/MaxBin2err.log'
    #    run(cmd)
    #    self.bins_as_fasta(f'{binner_dir}/output/bins')

    def run_main(self, binner_dir, options, threads, samples, binning_mode, target, coverage, coverage_dir, **kwargs):
        apc = getattr(cv, f'{coverage}Coverage').abundance_precomputed()
        makedirs(join(binner_dir, f'output/bins'), exist_ok=True)
        contigs = join(binner_dir, 'input', 'asm.fasta')
        coverage_dir = coverage_dir if apc else join(binner_dir, 'input')
        covs = [join(coverage_dir, f'target_{target}_coverage.tsv')] if (binning_mode == 'single' and apc) else [join(coverage_dir, f'{s}_coverage.tsv') for s in samples] 
        covs = ' '.join(f'-abund{i} {fl}' for i, fl in enumerate(covs,start=1))
        covs = covs.replace('-abund1', '-abund')
        cmd = f'{self.exec} -thread {threads} {options} -contig {contigs} {covs} -out {binner_dir}/output/bins/bin 2> {binner_dir}/output/MaxBin2err.log'
        run(cmd)
        self.bins_as_fasta(f'{binner_dir}/output/bins')

    #def run_prep(self, assembly_dir, map_dir, samples, binner_dir, **kwargs):
    #    target_sample = samples[0]
    #    makedirs(join(binner_dir, 'input'), exist_ok=True)
    #    # first get the contigs that will be binned (asm.fasta)
    #    softlink_assembly_to_taskdir(assembly_dir, target_sample, binner_dir)
    #    for s in samples:
    #        bam = join(map_dir, f'{target_sample}_{s}.bam')
    #        idxstats = join(binner_dir, 'input', f'{s}.idxstats')
    #        # compute idx stats
    #        compute_idxstats(bam, idxstats)
    #        counts = pd.read_csv(idxstats, sep='\t', header=None, index_col=0)[2]
    #        counts = counts[:-1]  # Drop a trailing '*' from the idxstats file
    #        counts.to_csv(idxstats.replace('.idxstats','.counts'), sep='\t', header=None)
    
    def run_prep(self, assembly_dir, coverage_dir, samples, binner_dir, target, coverage, **kwargs):
        apc = getattr(cv, f'{coverage}Coverage').abundance_precomputed()
        makedirs(join(binner_dir, 'input'), exist_ok=True)
        softlink_assembly_to_taskdir(assembly_dir, target, binner_dir)
        if apc: # nothing to do
            return
        bams = ' '.join([join(coverage_dir, f'{e}.bam') for e in samples])
        cov = join(binner_dir, 'input','coverage_mat_jgi.tsv')
        cv.jgi_summarize(bams, cov)
        cv.jgi_to_maxbin2(dirname(cov),cov)

    def bins_as_fasta(self, bin_dir):
        bins = glob.glob(join(bin_dir, '*.fasta'))
        return bins

    def prep_done(self, binner_dir, samples, coverage, **kwargs):
        fls = [join(binner_dir, 'input', f'{s}_coverage.tsv') for s in samples]
        return all(not_empty(fl) for fl in fls)

    def main_done(self, binner_dir, **kwargs):
        bin_dir = join(binner_dir, 'output/bins')
        return bin_done(bin_dir)

    
class METABAT2Binner(Binner):
    name='METABAT2'
    exec='metabat2'

    def run_main(self,binner_dir, options, threads, coverage_dir, binning_mode, coverage, **kwargs):
        apc = getattr(cv, f'{coverage}Coverage').abundance_precomputed()
        coverage_dir = coverage_dir if apc else join(binner_dir, 'input')
        makedirs(join(binner_dir, f'output/bins'), exist_ok=True)
        contigs = join(binner_dir, 'input', 'asm.fasta')
        cov = join(coverage_dir, ('target_' if (binning_mode == 'single' and apc) else '') + 'coverage_mat_jgi.tsv')
        cmd = f'{self.exec} --numThreads {threads} {options} --inFile {contigs}  --abdFile {cov} --outFile {binner_dir}/output/bins/bin'
        run(cmd)
        self.bins_as_fasta(f'{binner_dir}/output/bins')

    def run_prep(self, assembly_dir, samples, binner_dir, target, coverage, **kwargs):
        apc = getattr(cv, f'{coverage}Coverage').abundance_precomputed()
        makedirs(join(binner_dir, 'input'), exist_ok=True)
        softlink_assembly_to_taskdir(assembly_dir, target, binner_dir)
        if apc: # nothing to do
            return
        bams = ' '.join([join(assembly_dir, f'{e}.bam') for e in samples])
        cv.jgi_summarize(bams, join(binner_dir, 'input', 'coverage_mat_jgi.tsv'))

    def bins_as_fasta(self, bin_dir):
        bins = glob.glob(join(bin_dir, '*.fa'))
        if len(bins) != 0:
            for b in bins:
                rename(b, b.replace('.fa', '.fasta'))
        bins = glob.glob(join(bin_dir, '*.fasta'))
        return bins

    def prep_done(self, binner_dir, coverage, **kwargs):
        return not_empty(join(binner_dir, 'input', 'coverage_mat_jgi.tsv'))

    def main_done(self, binner_dir, **kwargs):
        bin_dir = join(binner_dir, 'output/bins')
        return bin_done(bin_dir)


class VAMBBinner(Binner):
    name='VAMB'
    exec='/insomnia001/depts/pmg/users/ic2465/miniforge3/envs/vamb/bin/vamb'

    def run_main(self, binner_dir, options, threads, coverage_dir, binning_mode, coverage, **kwargs):
        contigs = join(binner_dir, 'input', 'asm.fasta')
        makedirs(join(binner_dir, 'output'),exist_ok=True)
        bin_dir = join(binner_dir, 'output/bins')
        apc = getattr(cv, f'{coverage}Coverage').abundance_precomputed()
        if apc:
            pref = 'target_' if binning_mode == 'single' else ''
            abund = f'--abundance_tsv {join(coverage_dir, pref+'coverage_mat.tsv')}'
        else:
            abund = f'--bamdir {join(binner_dir, 'input')}'
        cmd = f'{self.exec} bin default -p {threads} {options} --outdir {bin_dir} --fasta {contigs} {abund} --minfasta 100000 -o "" 2> {binner_dir}/VAMB.errlog'
        run(cmd)
        self.bins_as_fasta(bin_dir)

    def run_prep(self, assembly_dir, coverage_dir, samples, binner_dir, target, coverage, **kwargs):
        makedirs(join(binner_dir, 'input'), exist_ok=True)
        softlink_assembly_to_taskdir(assembly_dir, target, binner_dir)
        apc = getattr(cv, f'{coverage}Coverage').abundance_precomputed()
        if apc: # nothing to do
            return
        for s in samples:
            bam = join(coverage_dir, f'{s}.bam')
            bamln = join(binner_dir, 'input', f'{s}.bam')
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

    def prep_done(self, binner_dir, samples, coverage, **kwargs):
        input_dir = join(binner_dir, 'input')
        return not_empty(f'{input_dir}/asm.fasta') and all(not_empty(f'{input_dir}/{s}.bam') for s in samples)

    def main_done(self, binner_dir, **kwargs):
        bin_dir = join(binner_dir, 'output/bins')
        return bin_done(bin_dir)
    
class SemiBin2Binner(Binner):
    name='SemiBin2'
    exec = 'SemiBin2'

    def run_main(self, binner_dir, options, samples, threads, coverage_dir, binning_mode, coverage, **kwargs):
        apc = getattr(cv, f'{coverage}Coverage').abundance_precomputed()
        bin_dir = join(binner_dir, 'output/bins')
        # single and multi-sample binning differ
        contigs = join(binner_dir, 'input', 'asm.fasta')
        #if len(samples) == 1:
        #    bam = join(binner_dir, 'input', f'{target_sample}_{target_sample}.bam')
        #    run(
        #        f'{self.exec} single_easy_bin --threads {threads} {options} --environment human_gut -i {contigs} -b {bam} ' + 
        #        f'-o {bin_dir} --compression none 2> {binner_dir}/SemiBin2.errlog'
        #    )
        #else:
        if apc: 
            pref = 'target_' if binning_mode == 'single' else ''
            if coverage == 'aemb':
                abund = '-a ' + ' '.join(f'{coverage_dir}/sb2_{s}_coverage.tsv' for s in 2*samples)
            else:
                abund = f'--depth-metabat2 {join(coverage_dir, f'{pref}coverage_mat_jgi.tsv')} --environment human_gut'
        else:
            abund = '-b ' + ' '.join(f'{binner_dir}/input/{s}.bam' for s in samples)
        cmd = f'{self.exec} single_easy_bin --threads {threads} {options} -i {contigs} {abund} ' + \
        f'-o {bin_dir} --compression none 2> {binner_dir}/SemiBin2.errlog'
        run(cmd)
        self.bins_as_fasta(bin_dir)

    def run_prep(self, assembly_dir, coverage_dir, samples, target, binner_dir, coverage, **kwargs):
        makedirs(join(binner_dir, 'input'), exist_ok=True)
        softlink_assembly_to_taskdir(assembly_dir, target, binner_dir)
        apc = getattr(cv, f'{coverage}Coverage').abundance_precomputed()
        if apc: # nothing to do
            return
        # single and multi-sample prep differ
        for sample in samples:
            bam = join(coverage_dir, f'{sample}.bam')
            bamln = join(binner_dir, 'input', f'{sample}.bam')
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

    def prep_done(self, binner_dir, samples, coverage, **kwargs):
        input_dir = join(binner_dir, 'input')
        return not_empty(f'{input_dir}/asm.fasta') and all(not_empty(f'{input_dir}/{s}.bam') for s in samples)

    def main_done(self, binner_dir, **kwargs):
        bin_dir = join(binner_dir, 'output/bins')
        return bin_done(bin_dir)
    

class CONCOCTBinner(Binner):
    # Fits directly into current framework include
    name='CONCOCT'
    exec='/insomnia001/depts/pmg/users/ic2465/miniforge3/envs/concoct/bin/'

    def run_main(self, binner_dir, options, threads, **kwargs):
        # run concoct main clustering
        input = join(binner_dir, 'input')
        bin_dir = join(binner_dir, 'output/bins')
        makedirs(bin_dir, exist_ok=True)
        run(
            f'{self.exec}/concoct --threads {threads} {options} --composition_file {input}/contig_10K.fa --coverage_file {input}/coverage_table.tsv -b {bin_dir} 2> {binner_dir}/output/cncterr.log'
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

    def run_prep(self, assembly_dir, coverage_dir, target, samples, binner_dir, **kwargs):
        makedirs(join(binner_dir, 'input'), exist_ok=True)
        # first get the contigs that will be binned
        asm = softlink_assembly_to_taskdir(assembly_dir, target, binner_dir)
        # cut up the fasta into chunks
        run(
            f'{self.exec}/cut_up_fasta.py {asm} -c 10000 -o 0 --merge_last -b {binner_dir}/input/contig_10K.bed ' +
            f'> {binner_dir}/input/contig_10K.fa'
        )
        # make the coverage table
        bams = ' '.join([join(coverage_dir, f'{target}_{s}.bam') for s in samples])
        run(
            f'{self.exec}/concoct_coverage_table.py {binner_dir}/input/contig_10K.bed {bams} ' +
            f'> {binner_dir}/input/coverage_table.tsv'
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

    def prep_done(self, binner_dir, **kwargs):
        input_dir = join(binner_dir, 'input')
        return not_empty(f'{input_dir}/contig_10K.bed') and not_empty(f'{input_dir}/contig_10K.fa') and not_empty(f'{input_dir}/coverage_table.tsv')

    def main_done(self, binner_dir, **kwargs):
        bin_dir = join(binner_dir, 'output/bins')
        return bin_done(bin_dir)


class COMEBinBinner:
    name='COMEBin'
    exec='comebin'
    size=1000
    def run_main(self, binner_dir, options, threads, **kwargs):
        output = join(binner_dir, 'output')
        input = join(binner_dir, 'input')
        rm_dir(join(output, 'bins', 'comebin_res'))
        run(
            f'mamba run -n {self.exec} python -c "import torch; print(torch.cuda.is_available(),flush=True)"'
        )
        run(
            f'mamba run -n {self.exec} run_comebin.sh -t {threads} {options} -a {input}/asm_{self.size}.fa ' + 
            f'-p {input} -o {output}/bins -s main 2> {output}/COMEBin.errlog'
        )
        self.bins_as_fasta(f'{output}/bins')

    def run_prep(self, assembly_dir, coverage_dir, samples, binner_dir, options, threads, **kwargs):
        input_dir = join(binner_dir, 'input')
        output_dir = join(binner_dir, 'output')
        makedirs(input_dir,exist_ok=True)
        makedirs(output_dir,exist_ok=True)
        makedirs(f'{output_dir}/bins', exist_ok=True)
        target_sample = samples[0]
        asm =  softlink_assembly_to_taskdir(assembly_dir, target_sample, binner_dir)
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
                f'samtools view -@ {threads} -b -L {input_dir}/contig_names{self.size}.bed {coverage_dir}/{target_sample}_{s}.bam ' + 
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

    def prep_done(self, binner_dir, **kwargs):
        input_dir = join(binner_dir, 'output')
        means = glob.glob(f'{input_dir}/bins/data_augmentation/*_mean.tsv')
        vars = glob.glob(f'{input_dir}/bins/data_augmentation/*_var.tsv')
        return len(means) == len(vars) and all(not_empty(e) for e in means+vars) and len(means) > 0

    def main_done(self, binner_dir, **kwargs):
        bin_dir = join(binner_dir, 'output/bins')
        return bin_done(bin_dir)

class Refiner(ABC):
    name: str
    exec: str

    @abstractmethod
    def run_main(self, binner_dir, options, **kwargs):
        ...

    @abstractmethod
    def main_done(self, binner_dir, **kwargs) -> bool:
        ...

    @abstractmethod
    def bins_as_fasta(self, bin_dir) -> list:
        ...

class DAS_ToolRefiner(Refiner):
    name = 'DAS_Tool'
    exec = 'DAS_Tool'
    def run_main(self, binner_dir, options, pipelines, ppln_bin_out,assembly_dir, target_sample, threads=8, **kwargs):
        summaries = list()
        bin_dir = join(binner_dir, f'output/bins')
        makedirs(bin_dir, exist_ok=True)
        input_dir = join(binner_dir, 'input') 
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
        contigs = softlink_assembly_to_taskdir(assembly_dir, target_sample, binner_dir)
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
    def main_done(binner_dir, **kwargs):
        bin_dir = join(binner_dir, 'output/bins')
        return bin_done(bin_dir)
