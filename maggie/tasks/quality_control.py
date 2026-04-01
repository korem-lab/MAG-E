from abc import ABC, abstractmethod
import pandas as pd
from os import makedirs, rename
from os.path import join, dirname
from subprocess import run
import glob

class QCTool(ABC):
    name: str
    exec: str
    qc_measures: list

    @abstractmethod
    def run(self, task_out_dir, threads, **kwargs):
        ...

    @abstractmethod 
    def to_qctable(self, task_out_dir, **kwargs):
        ...

class CheckM2QCTool(QCTool):
    name = 'CheckM2'
    exec = 'checkm2'
    qc_measures = [f'{name}_completeness', f'{name}_contamination']
    database_path = 'uniref100.KO.1.dmnd'

    def run(self, task_out_dir, threads, **kwargs):
        bin_dir = join(task_out_dir, 'output/bins')
        task_out_dir = join(task_out_dir, f'output/{self.name}')
        makedirs(task_out_dir, exist_ok=True)
        err_log = join(dirname(task_out_dir), 'checkM2err.log')
        cmd = f'{self.exec} predict -t {threads} --force -x fasta ' \
        f'--input {bin_dir} --output-directory {task_out_dir} --database_path {self.database_path} 2> {err_log}'
        run(cmd, shell=True)

    def to_qctable(self, task_out_dir, **kwargs):
        result_dir = join(task_out_dir, f'output/{self.name}')
        res = pd.read_csv(join(result_dir, 'quality_report.tsv'),sep='\t')
        res.rename({'Name':'bin', 'Completeness':'completeness','Contamination':'contamination'},axis=1,inplace=True)
        res= res.add_prefix(f'{self.name}_')
        res.rename({f'{self.name}_bin':'bin'},axis=1,inplace=True)
        return res
    
class GUNCQCTool(QCTool):
    name = 'GUNC'
    exec = 'gunc'
    database_path ='gunc_db_progenomes2.1.dmnd'
    qc_measures = [f'{name}_isPass']

    def run(self, task_out_dir, threads, **kwargs):
        bin_dir = join(task_out_dir, 'output/bins')
        task_out_dir = join(task_out_dir, f'output/{self.name}')
        makedirs(task_out_dir, exist_ok=True)
        run(
            f'{self.name} run -t {threads} --input_dir {bin_dir} --file_suffix fasta ' +
            f'-r {self.database_path} --out_dir {task_out_dir} 2> {task_out_dir}/GUNCerr.log', shell=True
        )

    def to_qctable(self, task_out_dir, **kwargs):
        results_dir = join(task_out_dir, 'output/GUNC')
        res = pd.read_csv(join(results_dir,'GUNC.progenomes_2.1.maxCSS_level.tsv'),sep='\t')
        res.rename({'genome':'bin','pass.GUNC':'isPass'},axis=1,inplace=True)
        res['bin'] = res['bin'].apply(lambda x: x[:-1] if x.endswith('.') else x)
        res = res.add_prefix(f'{self.name}_')
        res.rename({f'{self.name}_bin':'bin'},axis=1,inplace=True)
        return res