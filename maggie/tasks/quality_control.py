from abc import ABC, abstractmethod
import pandas as pd
from os import makedirs, rename
from os.path import join, dirname
from ..utils import not_empty, run
import glob

class QCTool(ABC):
    name: str
    exec: str
    qc_measures: list

    def has_bins(self, binner_dir):
        bins = glob.glob(join(binner_dir, 'output/bins/*.fasta'))
        return len(bins)>0

    @abstractmethod
    def run_main(self, binner_dir, threads, **kwargs):
        ...

    @abstractmethod 
    def cleanup(self, binner_dir, **kwargs):
        ...

    @abstractmethod 
    def main_done(self, binner_dir, **kwargs):
        ...

class CheckM2QCTool(QCTool):
    name = 'CheckM2'
    exec = '/insomnia001/depts/pmg/KoremLab/miniforge/envs/checkm2/bin/checkm2'
    qc_measures = [f'{name}_completeness', f'{name}_contamination']
    database_path = '/insomnia001/depts/pmg/KoremLab/Databases/checkm2/uniref100.KO.1.dmnd'

    def run_main(self, binner_dir, threads, **kwargs):
        bin_dir = join(binner_dir, 'output/bins')
        binner_dir = join(binner_dir, f'output/{self.name}')
        makedirs(binner_dir, exist_ok=True)
        err_log = join(dirname(binner_dir), 'checkM2err.log')
        cmd = f'{self.exec} predict -t {threads} --force -x fasta ' \
        f'--input {bin_dir} --output-directory {binner_dir} --database_path {self.database_path} 2> {err_log}'
        run(cmd)

    def cleanup(self, binner_dir, **kwargs):
        result_dir = join(binner_dir, f'output/{self.name}')
        res = pd.read_csv(join(result_dir, 'quality_report.tsv'),sep='\t')
        res.rename({'Name':'bin', 'Completeness':'completeness','Contamination':'contamination'},axis=1,inplace=True)
        res= res.add_prefix(f'{self.name}_')
        res.rename({f'{self.name}_bin':'bin'},axis=1,inplace=True)
        return res

    def prep_done(self, **kwargs):
        return True

    def main_done(self, binner_dir, **kwargs):
        results_dir = join(binner_dir, f'output/{self.name}')
        return not_empty(join(results_dir, 'quality_report.tsv'))
    
class GUNCQCTool(QCTool):
    name = 'GUNC'
    exec = 'gunc'
    database_path ='/insomnia001/depts/pmg/users/ab4966/gunc/gunc_db/gunc_db_progenomes2.1.dmnd'
    qc_measures = [f'{name}_isPass']

    def run_main(self, binner_dir, threads, **kwargs):
        bin_dir = join(binner_dir, 'output/bins')
        binner_dir = join(binner_dir, f'output/{self.name}')
        makedirs(binner_dir, exist_ok=True)
        cmd = f'{self.exec} run -t {threads} --input_dir {bin_dir} --file_suffix fasta -r {self.database_path} --out_dir {binner_dir} 2> {binner_dir}/GUNCerr.log'
        run(cmd)

    def cleanup(self, binner_dir, **kwargs):
        results_dir = join(binner_dir, f'output/{self.name}')
        res = pd.read_csv(join(results_dir,'GUNC.progenomes_2.1.maxCSS_level.tsv'),sep='\t')
        res.rename({'genome':'bin','pass.GUNC':'isPass'},axis=1,inplace=True)
        res['bin'] = res['bin'].apply(lambda x: x[:-1] if x.endswith('.') else x)
        res = res.add_prefix(f'{self.name}_')
        res.rename({f'{self.name}_bin':'bin'},axis=1,inplace=True)
        return res

    def prep_done(self, **kwargs):
        return True

    def main_done(self, binner_dir, **kwargs):
        results_dir = join(binner_dir, f'output/{self.name}')
        return not_empty(join(results_dir,'GUNC.progenomes_2.1.maxCSS_level.tsv'))