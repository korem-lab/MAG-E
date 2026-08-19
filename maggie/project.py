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
from .tasks import assembly as ab, binning as bn, quality_control as qc, task_utils as tu
from . import ground_truth as gt 
from .utils import parse_manifest, run_R_script, parse_contig_properties, print_and_return, soft_link
from .plotting import plot_pipeline_performance_model, plot_unlabelled_version


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
    def load(cls, project_path: Path, verify=True) -> "Project":
        config = Config.from_yaml(join(project_path, "config.yaml"))
        if verify:
            config.verify_config()
        return cls(project_path, config)

    def build_core_directories(self):
        # Verify the configuration 
        if self.config.simulate == 'yes':
            os.makedirs(self.config.ecosystem_db, exist_ok=True)
        else:
            # Soft link the reads
            fqs = glob.glob(join(self.config.samples_dir, '*.fastq.gz'))
            for srcfq in fqs:
                dstfq = join(self.config.simulation_dir, os.path.basename(srcfq))
                soft_link(srcfq, dstfq)
        os.makedirs(self.config.simulation_dir, exist_ok=True)
        os.makedirs(self.config.bintask_dir, exist_ok=True)
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
    
    def construct_tasks(self, write_to_disk=False):
        """
        Construct the manifest, which is a record of all MAG generation tasks.
        """

        # make the manifest
        manifest = tu.make_manifest(
            self.config.assemblers, self.config.mappers, self.config.binners, self.config.binning_modes, self.config.refiners, self.config.qctools,
            self.config.simulation_dir, self.config.assembly_cache, self.config.mapping_cache,
            self.config.bintask_dir, self.config.project_base
        )

        # make the core directories for MAG generation. 
        os.makedirs(self.config.assembly_cache, exist_ok=True)
        os.makedirs(self.config.bintask_dir, exist_ok=True)
        manifest.asm_dir.apply(lambda x: os.makedirs(x, exist_ok=True) if not pd.isna(x) else None)
        manifest.task_out_dir.apply(lambda x: os.makedirs(x, exist_ok=True))
        manifest.map_dir.apply(lambda x: os.makedirs(x,exist_ok=True) if not pd.isna(x) else None)

        if write_to_disk:
            manifest.to_csv(join(self.config.project_base, 'manifest.csv'), index=None)
        return manifest

    def get_tasks(
            self, manifest, assembler, aopt, mapper, mopt, binner, bopt, binning_mode, refiner, ropt, rpipe, target
        ):
        flt = lambda x,y: manifest[x] == y if y is not None else pd.Series(True, index=manifest.index)
        rpipe = tuple(tuple(e.split(',')) for e in rpipe.split(':'))
        return manifest.loc[
            flt('target_sample', target) & flt('assembler', assembler) & flt('assembler_options', aopt) &
            flt('mapper', mapper) & flt('mapper_options', mopt) & flt('binner', binner) & 
            flt('binner_options', bopt) & flt('binning_mode', binning_mode) &
            flt('refiner', refiner) & flt('refiner_opt', ropt) & flt('pipelines', rpipe)
        ]

    def get_task_directory(
            self, manifest, assembler, aopt, mapper, mopt, binner, bopt, binning_mode, refiner, ropt, rpipe, target
        ):
        tsks = self.get_tasks(
            manifest, assembler, aopt, mapper, mopt, binner, bopt, binning_mode, refiner, ropt, rpipe, target
        )
        if assembler and mapper and binner and binning_mode and target:
            assert (binner and not refiner) or (refiner and not binner), Exception('Pick binner xor refiner.')
            tsks[['task_out_dir']].drop_duplicates()
        elif assembler and mapper and target:
            tsks[['map_dir']].drop_duplicates()
        elif assembler and target:
            tsks[['asm_dir']].drop_duplicates()
        else:
            Exception("Invalid query.")
        if tsks.empty:
            raise Exception("No tasks fit the query.")
        elif len(tsks) > 1:
            raise Exception("More than task fits the query") 
        return tsks.iloc[0,0]

    def query(self, 
            type, target, assembler, aopt, binner, bopt, 
            binning_mode, mapper, mopt, refiner, ropt, rpipe, simulations,
            genomes, evaluations, of, within
        ):
        """
        Queries the maggie project for information
        """
        manifest = parse_manifest(self.config.manifest)
        assert type in ['dir', 'list'], Exception('Only "dir" and "list" are valid queries.')

        if type == 'dir':
            if simulations:
                return print_and_return(self.config.simulation_dir)
            if genomes:
                return print_and_return(self.config.genomes_dir)
            if evaluations:
                return print_and_return(self.config.evaluation_dir)
            task_dir = self.get_task_directory(
                manifest, assembler, aopt, mapper, mopt, binner, bopt, binning_mode, refiner, ropt, rpipe, target
            )
            return print_and_return(task_dir)
        elif type == 'list':
            assert of, Exception("--of must be specified if the query is a list.")
            assert of in ['binners', 'assemblers', 'mappers', 'modes', 'qctools', 'targets', 'options', 'samples'], Exception('Invalid --of argument.')
            if of in ['options', 'samples']:
                assert within, Exception("Asking for a list of options or samples requires within to be specified.")
            if of not in ['options', 'samples']:
                return print_and_return(self.config.get_list_of(of))
            elif of == 'options':
                return print_and_return(self.config.get_options_within(within))
            elif of == 'samples' and target:
                return print_and_return(self.config.get_samples_within_mode(within, target))
            else:
                raise Exception('Invalid list query.')
        else:
            raise Exception('Only "dir" and "list" are valid queries.')

    def run_task(
        self, type, index, target, assembler, aopt, binner, bopt, binning_mode, mapper, mopt, refiner, ropt, rpipe, threads
    ):
        """
        Generic interface to launch tasks MAG-E tasks from.
        """
        manifest = parse_manifest(self.config.manifest)
        if index:
            tsk = manifest.iloc[index,:]
        else:
            tsk = self.get_tasks(
                manifest, assembler, aopt, mapper, mopt, binner, bopt, binning_mode, refiner, ropt, rpipe, target
            )
        tsk = tsk if len(tsk) == 1 else tsk.loc[0,:]
        if type == 'assembly':
            tu.run_assembly(tsk)
        elif type == 'mapping':
            tu.run_mapping(tsk)
        elif type == 'binning':
            tu.run_binning(tsk)
        elif type == 'qc':
            tu.run_quality_control(tsk)

    def run_assembly(self, threads):
        """
        Runs each assembly task in the manifest. 
        """

        manifest = parse_manifest(self.config.manifest)
        tu.run_assemblies(manifest, threads)

    def run_mapping(self, threads):
        """
        Runs each assembly task in the manifest. 
        """

        manifest = parse_manifest(self.config.manifest)
        tu.run_mapping(manifest, threads)

    def construct_ground_truth(self, min_contig_len=100, min_pident=99, min_prop=99, max_prop=101):
        """
        Constructs the ground truth for each assembly. 
        """
        manifest = parse_manifest(self.config.manifest)
        gt.construct_ground_truth(
            manifest, self.config.ecosystem_db_metadata,
            min_contig_len, min_pident, min_prop, max_prop
        )
    
    def run_binning(self, run_only, force_prep, force_bin, threads=8):
        manifest = parse_manifest(self.config.manifest)
        tu.run_bin_tasks(manifest, run_only, force_prep, force_bin, threads)

    def run_refine(self, run_only, force_refine, threads=8):
        manifest = parse_manifest(self.config.manifest)
        tu.run_refine_tasks(manifest, run_only, force_refine, threads)
    
    def run_quality_control(self, threads=8):
        manifest = parse_manifest(self.config.manifest)
        tu.run_quality_control(manifest, threads=threads, force=True)

    def calc_per_genome_metrics(self):
        # Construct the binning table for each binning tasks
        manifest = parse_manifest(self.config.manifest)
        tu.construct_binning_tables(manifest)

        # Build the reports. There is a report for each assembler, binner pair. 
        # A report contains the relevant information (genome, abundance, assigned bin, bin quality control)
        # for all (assembler, binner)-combo bin tasks. These reports are used to calculate the genome metrics.
        os.makedirs(self.config.evaluation_dir, exist_ok=True)
        bnab = manifest[['binner', 'assembler', 'is_refiner']].dropna().drop_duplicates().values
        rfab = manifest[['refiner', 'assembler', 'is_refiner']].dropna().drop_duplicates().values
        combos = np.vstack([bnab, rfab])
        for i in range(combos.shape[0]):
            bn, ab, is_refiner = combos[i]
            if is_refiner:
                _mnfst = manifest.loc[(manifest.refiner == bn) & (manifest.assembler == ab)]
            else:
                _mnfst = manifest.loc[(manifest.binner == bn) & (manifest.assembler == ab)]
            rprt = ev.construct_report(_mnfst, add_cp=False, read_counts=self.config.read_counts)
            rprt.to_parquet(join(self.config.evaluation_dir, f'{ab}_{bn}_REPORT.parquet'))

        # Compute the per-genome metrics
        for i in range(combos.shape[0]):
            bn, ab, _ = combos[i]
            rprt = pd.read_parquet(join(self.config.evaluation_dir, f'{ab}_{bn}_REPORT.parquet'))
            gm = ev.construct_genome_metrics(rprt)
            # these got dropped when constructing the genome metrics, all were doing here is adding them back
            gm = ev.add_report_data(rprt, gm, manifest.copy())
            gm.to_parquet(join(self.config.evaluation_dir, f'{ab}_{bn}.genome_metrics.parquet'))

    def calc_contig_level_metrics(self, precision, recall, drop_raw_property=True):
        manifest = parse_manifest(self.config.manifest)
        # Add ground truth genome origin to the contig properties
        asm_tasks = manifest.loc[['assembler', 'target_sample', 'assembly_options', 'asm_dir']].drop_duplicates()
        asm_tasks = asm_tasks[asm_tasks.assembly_options == 'default']
        cp_df_dataset = list()
        for i in range(len(asm_tasks)):
            t = asm_tasks.iloc[i,:]
            cp_df = parse_contig_properties(join(t.asm_dir, f'{t.target_sample}_contig_properties.csv'))
            gt_df = pd.read_csv(join(t.asm_dir, f'{t.target_sample}_gt_table.csv'))
            cp_df = pd.merge(cp_df, gt_df[['genome', 'key']], on='key')
            cp_df_dataset.append(cp_df)
        cp_df_dataset = pd.concat(cp_df_dataset)
        cp_df_dataset.reset_index(drop=True, inplace=True)

        # collect the genome metrics over the default tasks
        dflt_bin_tasks = manifest[
            (manifest.binner_options == 'default') & (manifest.assembly_options == 'default')
        ]
        gnm_metrics = ev.parse_genome_measurement_files(
            join(self.config.evaluation_dir, f'*.genome_metrics.parquet')
        )
        gnm_metrics = gnm_metrics.loc[gnm_metrics.task_name.isin(dflt_bin_tasks.task_name)]

        # get the recoverable set and filter to just those 
        _, rcvgnms = ev.recoverable_genome_set(gnm_metrics, recall, precision)
        cp_df_dataset = pd.merge(cp_df_dataset, rcvgnms, on=['sample', 'genome'])

        # iterate over continuous propertis and make a percentile form
        for prop in [e for e in cp_df_dataset.columns if e.startswith('prop_c')]:
            cp_df_dataset = pd.merge(
                cp_df_dataset, ev.to_percentiles(cp_df_dataset[[prop, 'contig_length']]),
                left_index=True, right_index=True, how='left'
            )
            if drop_raw_property:
                cp_df_dataset.drop(prop, axis=1, inplace=True)
        cp_df_dataset.to_parquet(self.config.evaluation_dir, 'processed_contig_properties.parquet')

        # break the eval up by binner and assembler combos
        bnab = manifest[['binner', 'assembler', 'is_refiner']].dropna().drop_duplicates().values
        rfab = manifest[['refiner', 'assembler', 'is_refiner']].dropna().drop_duplicates().values
        combos = np.vstack([bnab, rfab])
        for i in range(combos.shape[0]):
            bn, ab, is_refiner = combos[i]
            if is_refiner:
                _mnfst = manifest.loc[(manifest.refiner == bn) & (manifest.assembler == ab)]
            else:
                _mnfst = manifest.loc[(manifest.binner == bn) & (manifest.assembler == ab)]
            rprt = ev.construct_report(_mnfst, add_cp=True, add_abundance=False)
            rprt.to_parquet(join(self.config.evaluation_dir, f'{ab}_{bn}_CPREPORT.parquet'))
            ev.construct_contig_property_metrics(
                join(self.config.evaluation_dir, f'{ab}_{bn}_CPREPORT.parquet'), rcvgnms, manifest
            )

    def evaluate_pipelines(self, precision, recall, plots=False):
        # Write the set of recoverable genome scores to a csv for analysis in R
        os.makedirs(join(self.config.evaluation_dir, 'LMM'), exist_ok=True)
        gnm_metrics = ev.parse_genome_measurement_files(
            glob.glob(join(self.config.evaluation_dir, '*.genome_metrics.parquet')), min_cov=0
        )
        scores, _ = ev.recoverable_genome_set(gnm_metrics, 0.2, 0.2)
        scores = scores[['binning_mode', 'binner', 'assembler', 'genome', 'is_refiner', 'sample', 'metric', 'value']]
        scores.to_csv(join(self.config.evaluation_dir, 'LMM/per_genome_metrics.csv'), index=None)

        # Run the LMM and write output to disk
        run_R_script('LMM_optimal_mag_pipeline', join(self.config.evaluation_dir, 'LMM'), 'true')
        ## plot data if needed
        #if plots:
        #    for by in ['binner', 'binning-mode', 'assembler']:
        #        for metric in ['fscore', 'precision', 'recall']:
        #            order = ['MaxBin2', 'VAMB', 'CONCOCT', 'METABAT2', 'SemiBin2', 'COMEBin']
        #            ax = plot_pipeline_performance_model(
        #                pd.read_csv(join(self.config.evaluation_dir, 'LMM', f'all_pipelines_{metric}.csv')), by=by, y_name=metric, order=order
        #            )
        #            plot_unlabelled_version(ax, join(utls.MANU_FIGS, 'Fig2', 'MAG_pipeline_LMM_fscore_binner'),ylower=-0.05, yupper=1.3,show=False)
