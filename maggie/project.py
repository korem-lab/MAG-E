import pandas as pd
import os
import sys
from os.path import join, exists
import glob
import numpy as np
from pathlib import Path
from itertools import product

from .config import Config
from . import database as db
from . import mirror_specs as ms
from . import simulation as sm
from . import evaluation as ev
from .tasks import assembly as ab, binning as bn, quality_control as qc, coverage as cv, task_utils as tu
from . import ground_truth as gt 
from .utils import run_R_script, parse_contig_properties, print_and_return, soft_link, parse_read_counts
from .plotting import plot_pipeline_performance_model, plot_unlabelled_version
from .manifest import Manifest
from .evaluation import construct_binning_table, construct_quality_control_tables, construct_genome_metrics


class Project:

    def __init__(self, path: Path, config: Config):
        self.path = path
        self.config = config

    @classmethod
    def create_empty(cls, name: str, path: Path, use_api: str):
        project_path = join(path ,name)
        os.makedirs(project_path, exist_ok=True)
        config = Config(project_name=name, project_base=str(project_path), use_api=use_api)
        config.to_yaml(project_path)

    @classmethod
    def create(cls, name: str, path: str, config: Config) -> "Project":
        project_path = join(path, name)
        project_path.mkdir(parents=True)
        return cls(project_path, config) 

    @classmethod
    def load(cls, project_path: Path) -> "Project":
        config = Config.from_yaml(join(project_path, "config.yaml"))
        config.verify_config()
        return cls(project_path, config)

    def build_core_directories(self):
        if self.config.simulate == 'yes':
            os.makedirs(self.config.ecosystem_db, exist_ok=True)
        else:
            # Soft link the reads
            os.makedirs(self.config.simulation_dir, exist_ok=True)
            fqs = glob.glob(join(self.config.samples_dir, '*.fastq.gz'))
            for srcfq in fqs:
                dstfq = join(self.config.simulation_dir, os.path.basename(srcfq))
                soft_link(srcfq, dstfq)
        os.makedirs(self.config.simulation_dir, exist_ok=True)
        os.makedirs(self.config.evaluation_dir, exist_ok=True)

    def make_database(self, threads=8, ani=0.98, c=200, write_to_disk=False):
        """
        Makes the MAG-E database from the input genomes. 
        
        Each real sample is matched against the genomes of this database. Close
        genome matches are used to construct the mirror specification of the sample.
        """
        # make the directory housing the MAG-E database
        os.makedirs(self.config.ecosystem_db, exist_ok=True)

        # the database table will be built up from the user-specified cluster assignments
        db_table = pd.read_csv(self.config.ecosystem_db_metadata)
        db.check_genome_files_exist(db_table.genome, self.config.genomes_dir)

        # Cluster the genomes at the strain level
        db.run_strain_clustering(db_table, self.config.genomes_dir, self.config.strain_dreps, threads, ani)

        # Add the strain cluster information to the database table
        db_table = db.build_database_table(db_table, self.config.genomes_dir, self.config.strain_dreps)

        # build syldb of species representatives
        reprs = db_table.FileLocation[db_table.isSpeciesRepr]
        db.construct_sylphdb(reprs, self.config.ecosystem_db, db_prefix='repr', t=threads, c=c, force=False)

        # build sylphdb of all genomes
        all_genomes = db_table.FileLocation
        db.construct_sylphdb(all_genomes, self.config.ecosystem_db, db_prefix='all', t=threads, c=c, force=False)

        if write_to_disk:
            db_table.to_csv(self.config.ecosystem_db_metadata, index=None)
        return db_table
    
    def make_mirrors(self, threads, c, seed):
        """
        Constructs the mirror specifications for each sample. 
        """

        # Make a Sylph sketch (paired.sylsp) of each sample.
        ms.construct_sylphsp(self.config.samples_dir, self.config.sylsp_dir, self.config.prefix1, t=threads,c=c)

        ## Sylph query and profile each sample.
        ms.sylph_profile(self.config.sylsp_dir, self.config.ecosystem_db, self.config.simulation_dir, 'repr', t=threads)
        ms.sylph_query(self.config.sylsp_dir, self.config.ecosystem_db, self.config.simulation_dir, 'all', t=threads)

        # Collect the Sylph results, and construct the mirror specifications.
        profiles = sorted(glob.glob(join(self.config.simulation_dir, '*_sylph_profile.tsv')))
        queries = sorted(glob.glob(join(self.config.simulation_dir, '*_sylph_query.tsv')))
        for profile, query in zip(profiles, queries):
            ms.construct_metagenomic_specification(self.config.simulation_dir, self.config.ecosystem_db_metadata, profile, query, seed)

    def simulate_mgx(self, threads, n_reads, seed, print_script=False):
        """
        Simulated metagenomic data. 
        """
        specs = sorted(glob.glob(join(self.config.simulation_dir, '*_metagenome_spec.csv')))
        for spec in specs:
            sm.run_InSilicoSeq(
                spec, self.config.simulation_dir, self.config.read_counts, 
                n_reads, threads, seed=seed, force=False, print_script=print_script
            )
        if print_script:
            return 
        # Final iss cleanup
        fls = glob.glob(f'{self.config.simulation_dir}*iss.tmp*')
        for fl in fls:
            os.remove(fl)
    
    def construct_tasks(self):
        """
        Construct the manifest, which is a record of all MAG generation tasks.
        """

        # make the manifest
        manifest = Manifest(self.config.project_base)
        manifest.sync(self.config)
        # make the core directories for MAG generation. 
        os.makedirs(self.config.assembly_cache, exist_ok=True)
        os.makedirs(self.config.bintask_dir, exist_ok=True)
        manifest.m.assembler_dir.apply(lambda x: os.makedirs(x, exist_ok=True) if not pd.isna(x) else None)
        manifest.m.binner_dir.apply(lambda x: os.makedirs(x, exist_ok=True))
        manifest.m.coverage_dir.apply(lambda x: os.makedirs(x,exist_ok=True) if not pd.isna(x) else None)
        return manifest

    def query(self, 
            type, target, of, within, simulations, genomes, evaluations, **kwargs
        ):
        """
        Queries the maggie project for information
        """
        manifest = Manifest(self.config.project_base)
        assert type in ['dir', 'list'], Exception('Only "dir" and "list" are valid queries.')

        if type == 'dir':
            if simulations:
                return print_and_return(self.config.simulation_dir)
            if genomes:
                return print_and_return(self.config.genomes_dir)
            if evaluations:
                return print_and_return(self.config.evaluation_dir)
            task_dir = manifest.get_task_directory(target=target, **kwargs)
            return print_and_return(task_dir)
        elif type == 'list':
            assert of, Exception("--of must be specified if the query is a list.")
            assert of in ['binners', 'assemblers', 'refiners', 'coverage', 'modes', 'qctools', 'targets', 'options', 'samples'], Exception('Invalid --of argument.')
            if of in ['options', 'samples']:
                assert within, Exception("Asking for a list of options or samples requires within to be specified.")
            if of not in ['options', 'samples']:
                return print_and_return(self.config.get_list_of(of))
            elif of == 'options':
                return print_and_return(self.config.get_options_within(within))
            elif of == 'binner_sets':
                return print_and_return(self.config.get_binner_sets_within(**kwargs))
            elif of == 'samples' and target:
                return print_and_return(self.config.get_samples_within_mode(within, target))
            else:
                raise Exception('Invalid list query.')
        else:
            raise Exception('Only "dir" and "list" are valid queries.')

    def report(self, type, stage, _tool, to_file=False):
        assert stage in ['prep', 'main', 'all']
        manifest = Manifest(self.config.project_base)
        records = list()
        for idx in manifest.m.index:
            tsk = manifest.m.loc[idx,:]
            # other than qctools each tsk has just one tool per type
            # but we need to account for qctools having multiple tools per task
            tools = {
                'assembler': [tsk.assembler],
                'coverage': [tsk.coverage],
                'binner': [tsk.binner],
                'qctool': tsk.qctools
            }[type]
            for tool in tools:
                if _tool is not None and tool != _tool:
                    continue
                o = getattr(ab, f'{tool}Assembler')() if type == 'assembler' \
                    else getattr(cv, f'{tool}Coverage')() if type == 'coverage' \
                    else getattr(bn, f'{tool}Binner')() if type == 'binner' \
                    else getattr(qc, f'{tool}QCTool')()
                done = True
                if stage == 'prep' or stage == 'all':
                    done &= o.prep_done(**tsk.to_dict())
                if stage == 'main' or stage == 'all':
                    done &= o.main_done(**tsk.to_dict())
                if not done:
                    path = tsk.assembler_dir if type == 'assembler' else tsk.coverage_dir if type == 'coverage' else tsk.binner_dir
                    opts = tsk.assembler_options if type == 'assembler' else tsk.coverage_options if type == 'coverage' else tsk.binner_options if type == 'binner' else 'default'
                    records.append((tool, opts, tsk.target, type, path))
        records = pd.DataFrame(records, columns = ['tool', 'options', 'target', 'type', 'path']).drop_duplicates()
        if to_file:
            records.to_csv(
                os.path.join(
                    self.config.project_base, 'task_report.csv'
                ), index=None
            )
        else:
            print(records.to_string(), flush=True)
            
    def run_task(
        self, target, assembler, aopt, coverage, copt, cov_sample, binner, bopt, binning_mode, binner_set, qctool, stage, force, threads, check
    ):
        """
        Generic interface to launch tasks MAG-E tasks from.
        """
        manifest = Manifest(self.config.project_base)
        tsk = manifest.get_tasks(
            assembler, aopt, coverage, copt, binner, bopt, 
            binning_mode, binner_set, target, cov_sample
        )
        if tsk.empty:
            raise Exception("No tasks fit the query.")
        if assembler and coverage and binner and binning_mode:
            if (len(tsk)!=1): raise Exception("Ambiguous. More than one task possible..")
            if qctool:
                assert qctool in tsk.qctools.iloc[0]
                return not tu.run_quality_control(tsk.iloc[0,:], qctool, threads, force, check)
            else:
                return not tu.run_binning(tsk.iloc[0,:], stage, threads, force, check)

        elif assembler and coverage and not binner and not qctool:
            assert(len(tsk[['assembler', 'assembler_options','coverage','coverage_options','target']].drop_duplicates()) == 1)
            if cov_sample:
                # We're calculating the coverage for a specified sample (cov_sample) against the target
                return not tu.run_coverage(tsk.iloc[0,:], stage, cov_sample, threads, force, check)
            else:
                # We're calculating coverage matrices for the target. Calculate the most expansive coverage matrix:
                # if "all" mode, use "all" task, else use "single" task.
                if any(tsk.binning_mode == 'all'):
                    return not tu.run_coverage(tsk.loc[tsk.binning_mode == 'all',:].iloc[0,:], stage, cov_sample, threads, force, check)
                elif any(tsk.binning_mode == 'single'):
                    return not tu.run_coverage(tsk.loc[tsk.binning_mode == 'single',:].iloc[0,:], stage, cov_sample, threads, force, check)
                else:
                    raise Exception("Must have either all or single binning modes specified, but none were found in manifest.")

        elif assembler and not coverage and not binner and not qctool:
            assert(len(tsk[['assembler', 'assembler_options','target']].drop_duplicates()) == 1)
            return not tu.run_assembly(tsk.iloc[0,:], threads, force, check)
        else:
            raise Exception("Invalid options.")

    def construct_ground_truth(
            self, target, assembler, assembler_options, 
            min_contig_len=100, min_pident=99, min_prop=99, max_prop=101, threads=8
        ):
        """
        Constructs the ground truth for each assembly. 
        """
        manifest = Manifest(self.config.project_base)
        if target is not None and assembler is not None:
            tsks = manifest.get_tasks(assembler, aopt=assembler_options, target=target)
            assert len(tsks[['target', 'simulation_dir', 'assembler_dir']].drop_duplicates()) == 1
        else:
            tsks = manifest.m
        gt.construct_ground_truth(
            tsks, self.config.ecosystem_db_metadata,
            min_contig_len, min_pident, min_prop, max_prop, threads
        )

    def calculate_metrics(self,
            level, target, assembler, assembler_options, coverage, coverage_options, 
            binner, binner_options, binning_mode, binner_set, task_hash
        ):
        assert level in ['genome-level', 'contig-level'], Exception('level must be "genome-level" or "contig-level".')
        manifest = Manifest(self.config.project_base)
        if task_hash:
            task = manifest.get_tasks(task_hash)
        else:
            task = manifest.get_tasks(assembler, assembler_options, coverage, coverage_options, binner, binner_options, binning_mode, binner_set, target)
        assert len(task) == 1, Exception("Invalid task specification: more than one task matched description.")
        task = task.iloc[0]

        if level == 'genome-level':
            bintbl = construct_binning_table(task)
            qctbls = construct_quality_control_tables(task)
            read_count = parse_read_counts(self.config.read_counts)
            read_count = read_count.loc[target].item()
            ecodb = pd.read_csv(self.config.ecosystem_db_metadata)
            spec = pd.read_csv(join(self.config.simulation_dir, f'{target}_metagenome_spec.csv'))
            genome_metrics = construct_genome_metrics(bintbl, qctbls, read_count, ecodb, spec)
            asm = self.get_name(task.assembler, task.assembler_summary_name)
            cov = self.get_name(task.coverage, task.coverage_summary_name)
            bin = self.get_name(task.binner, task.binner_summary_name, task.binner_set)
            genome_metrics = genome_metrics.assign(
                assembler=asm, coverage=cov, binner=bin, binning_mode=binning_mode, task_hash=task_hash
            )
            genome_metrics.to_parquet(
                join(task.binner_dir, 'genome_metrics.pqt')
            )

    def get_name(self, tool, sname, bs=None):
        bs='('+'-'.join(bs)+')' if bs else ''
        return tool if sname == 1 else f'{tool}({sname}){bs}'

    def run_evaluation(self, level):
        manifest = Manifest(self.config.project_base)
        if level == 'genome-level':
            # Concat all the genome metrics 
            if exists(join(self.config.evaluation_dir, 'genome_metrics.pqt')):
                mets = pd.read_parquet(join(self.config.evaluation_dir, 'genome_metrics.pqt'))
            else:
                mets = pd.concat(
                    [pd.read_parquet(join(d, 'genome_metrics.pqt')) for d in manifest.m.binner_dir]
                )
                mets.to_parquet(self.config.evaluation_dir, 'genome_metrics.pqt')

            #llmcalls
