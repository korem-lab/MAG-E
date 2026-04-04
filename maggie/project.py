import pandas as pd
import os
from os.path import join, exists
import glob
from pathlib import Path
from itertools import product

from .config import Config
from . import database as db
from . import mirror_specs as ms
from . import simulation as sm
from . import evaluation as ev
from .tasks import assembly as ab, binning as bn, quality_control as qc, task_utils as tu
from . import ground_truth as gt 
from .utils import manifest_get, parse_manifest, run_R_script
from .plotting import plot_pipeline_performance_model, plot_unlabelled_version


class Project:

    def __init__(self, path: Path, config: Config):
        self.path = path
        self.config = config

    @classmethod
    def create_empty(cls, name: str, path: Path, use_api: str):
        project_path = join(path ,name)
        project_path.mkdir(parents=True)
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
        os.makedirs(self.config.maggie_db_dir, exist_ok=True)
        os.makedirs(self.config.simulation_dir, exist_ok=True)
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
        os.makedirs(self.config.maggie_db_dir, exist_ok=True)

        # the database table will be built up from the user-specified cluster assignments
        db_table = pd.read_csv(self.config.cluster_assignments)
        db.check_genome_files_exist(db_table.genome, self.config.genomes_dir)

        # Cluster the genomes at the strain level
        db.run_strain_clustering(db_table, self.config.genomes_dir, self.config.strain_dreps, threads, ani)

        # Add the strain cluster information to the database table
        db_table = db.build_database_table(db_table, self.config.genomes_dir, self.config.strain_dreps)

        # build syldb of species representatives
        reprs = db_table.FileLocation[db_table.isSpeciesRepr]
        db.construct_sylphdb(reprs, self.config.maggie_db_dir, db_prefix='repr', t=threads, c=c, force=True)

        # build sylphdb of all genomes
        all_genomes = db_table.FileLocation
        db.construct_sylphdb(all_genomes, self.config.maggie_db_dir, db_prefix='all', t=threads, c=c, force=True)

        if write_to_disk:
            db_table.to_csv(self.config.maggie_db_md, index=None)
        return db_table
    
    def make_mirrors(self, threads, c, seed):
        """
        Constructs the mirror specifications for each sample. 
        """

        # Make a Sylph sketch (paired.sylsp) of each sample.
        ms.construct_sylphsp(self.config.samples_dir, self.config.sylsp_dir, self.config.prefix1, t=threads,c=c)

        ## Sylph query and profile each sample.
        ms.sylph_profile(self.config.sylsp_dir, self.config.maggie_db_dir, self.config.simulation_dir, 'repr', t=threads)
        ms.sylph_query(self.config.sylsp_dir, self.config.maggie_db_dir, self.config.simulation_dir, 'all', t=threads)

        # Collect the Sylph results, and construct the mirror specifications.
        profiles = sorted(glob.glob(join(self.config.simulation_dir, '*_sylph_profile.tsv')))
        queries = sorted(glob.glob(join(self.config.simulation_dir, '*_sylph_query.tsv')))
        for profile, query in zip(profiles, queries):
            ms.construct_metagenomic_specification(self.config.simulation_dir, self.config.maggie_db_md, profile, query, seed)

    def simulate_mgx(self, threads, n_reads, seed):
        """
        Simulated metagenomic data. 
        """
        specs = sorted(glob.glob(join(self.config.simulation_dir, '*_metagenome_spec.csv')))
        for spec in specs:
            sm.run_InSilicoSeq(
                spec, self.config.simulation_dir, self.config.read_counts, 
                n_reads, threads, seed=seed, force=False
            )
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
            self.config.assemblers, self.config.binners, self.config.binning_modes, self.config.refiners, self.config.qctools,
            self.config.simulation_dir, self.config.assembly_cache, self.config.bintask_dir, self.config.project_base
        )

        # make the core directories for MAG generation. 
        os.makedirs(self.config.assembly_cache, exist_ok=True)
        os.makedirs(self.config.bintask_dir, exist_ok=True)
        manifest.asm_dir.apply(lambda x: os.makedirs(x, exist_ok=True) if not pd.isna(x) else None)
        manifest.task_out_dir.apply(lambda x: os.makedirs(x, exist_ok=True))

        if write_to_disk:
            manifest.to_csv(join(self.config.project_base, 'manifest.csv'), index=None)
        return manifest
    
    def get_task_dir(self, request, sample, assembler, assembler_options, binner, binner_options, refiner, refiner_options, binning_mode, pipelines):
        """
        Prints the assembly or bin task director for a task given the options. 
        """

        manifest = parse_manifest(self.config.manifest)
        assert request in ['assembly', 'bin'], Exception(f'Request must be either "assembly" or "bin", not {request}.')
        assert (manifest.target_sample == sample).any(), Exception(f'Sample {sample} not in target samples')


        if request == 'assembly':
            assert (manifest.assembler == assembler).any(), Exception(f'Assembler {assembler} not specified in config.')
            assert (manifest.assembler_options == assembler_options).any(), Exception(f'Assembler options {assembler_options} not specified in config.')
            tsks = manifest.loc[
                (manifest['target_sample'] == sample) & (manifest['assembler'] == assembler) & (manifest['assembler_options'] == assembler_options)
            ][['asm_dir']].drop_duplicates()
        else:
            assert ((assembler and binner and binning_mode) or (refiner and pipelines)), Exception('Must specify an assembler, binner and binning_mode, or a refiner with pipelines')
            if binner and refiner is None:
                assert (manifest.assembler ==assembler).any(), Exception(f'Assembler {assembler} not specified in config.')
                assert (manifest.assembler_options == assembler_options).any(), Exception(f'Assembler options {assembler_options} not specified in config.')
                assert (manifest.binner == binner).any(), Exception(f'Binner {binner} not specified in config.')
                assert (manifest.binner_options == binner_options).any(), Exception(f'Binner options {binner_options} not specified in config.')
                assert (manifest.binning_mode == binning_mode).any(), Exception(f'Binning mode {binning_mode} not specified in config.')

                tsks = manifest.loc[
                    (manifest['target_sample'] == sample) & (manifest['assembler'] == assembler) & (manifest['assembler_options'] == assembler_options) &
                    (manifest['binner'] == binner) & (manifest['binner_options'] == binner_options) & (manifest['binning_mode'] == binning_mode)
                ][['task_out_dir']]
            elif binner is None and refiner: 
                assert (manifest.refiner == refiner).any(), Exception(f'Refiner {refiner} not specified in config.')
                assert (manifest.refiner_options == refiner_options).any(), Exception(f'Refiner options {refiner_options} not specified in config.')

                pipelines = pipelines.split(':')
                pipelines = tuple(tuple(e.split(',')) for e in pipelines)
                tsks = manifest.loc[
                    (manifest['target_sample'] == sample) & (manifest['refiner'] == refiner) & 
                    (manifest['refiner_options'] == refiner_options) & (manifest['pipelines'] == pipelines)
                ][['task_out_dir']]
            else:
                raise Exception("Invalid combination.")

        if tsks.empty:
            raise Exception("Invalid combination. There are no tasks.")
        elif len(tsks) > 1:
            raise Exception("Input combination is underspecified. More than one assembly tasks shares that combination.")
        else:
            print(tsks.iloc[0,0], flush=True, end='')
            return(tsks.iloc[0,0])

    def run_assembly(self, threads):
        """
        Runs each assembly task in the manifest. 
        """

        manifest = parse_manifest(self.config.manifest)
        tu.run_assemblies(manifest, threads)

    def construct_ground_truth(self, min_contig_len=100, min_pident=99, min_prop=99, max_prop=101):
        """
        Constructs the ground truth for each assembly. 
        """
        manifest = parse_manifest(self.config.manifest)
        gt.construct_ground_truth(
            manifest, self.config.maggie_db_md,
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
        tu.run_quality_control(manifest, threads)

    def calc_per_genome_metrics(self):
        # Construct the binning table for each binning tasks
        manifest = parse_manifest(self.config.manifest)
        tu.construct_binning_tables(manifest)

        # Build the reports. There is a report for each assembler, binner pair. 
        # A report contains the relevant information (genome, abundance, assigned bin, bin quality control)
        # for all (assembler, binner)-combo bin tasks. These reports are used to calculate the genome metrics.
        binners = manifest.binner.unique()
        assemblers = manifest.assemblers.unique()
        os.makedirs(self.config.evaluation_dir)
        for bn, ab in product(binners, assemblers):
            _mnfst = manifest.loc[(manifest.binner == bn) & (manifest.assembler == ab)]
            rprt = ev.construct_report(_mnfst, add_cp=False)
            rprt.to_parquet(join(self.config.evaluation_dir, f'{ab}_{bn}_REPORT.parquet'))

        # Compute the per-genome metrics
        for bn, ab in product(binners, assemblers):
            rprt = pd.read_parquet(join(self.config.evaluation_dir, f'{ab}_{bn}_REPORT.parquet'))
            gm = ev.construct_genome_metrics(rprt)
            # these got dropped when constructing the genome metrics, all were doing here is adding them back
            gm = ev.add_report_data(rprt, gm, manifest)
            gm.to_parquet(join(self.config.evaluation_dir, f'{ab}_{bn}.genome_metrics.parquet'))
    
    def evaluate_pipelines(self, precision, recall, plots=False):
        # Write the set of recoverable genome scores to a csv for analysis in R
        os.makedirs(join(self.config.evaluation_dir, 'LMM'))
        gmfiles = glob.glob(self.config.evaluation_dir, '*.genome_metrics.parquet')
        scores, _ = ev.recoverable_genome_set(gmfiles, recall, precision)
        scores = scores[['binning_mode', 'binner', 'assembler', 'genome', 'sample', 'metric', 'value']]
        scores.to_csv(join(self.config.evaluation_dir, 'LMM/per_genome_metrics.csv'), index=None)

        # Run the LMM and write output to disk
        run_R_script('LMM_optimal_mag_pipeline', join(self.config.evaluation_dir, 'LMM'), 'true')
        # plot data if needed
        if plots:
            for by in ['binner', 'binning-mode', 'assembler']:
                for metric in ['fscore', 'precision', 'recall']:
                    order = ['MaxBin2', 'VAMB', 'CONCOCT', 'METABAT2', 'SemiBin2', 'COMEBin']
                    ax = plot_pipeline_performance_model(
                        pd.read_csv(join(self.config.evaluation_dir, 'LMM', f'all_pipelines_{metric}.csv')), by=by, y_name=metric, order=order
                    )
                    plot_unlabelled_version(ax, join(utls.MANU_FIGS, 'Fig2', 'MAG_pipeline_LMM_fscore_binner'),ylower=-0.05, yupper=1.3,show=False)
