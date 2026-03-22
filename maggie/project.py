import pandas as pd
import os
from os.path import join, exists
import glob
from pathlib import Path

from .config import Config
from . import database as db
from . import mirror_specs as ms


class Project:

    def __init__(self, path: Path, config: Config):
        self.path = path
        self.config = config

    @classmethod
    def create_empty(cls, name: str, path: Path):
        project_path = path / name
        project_path.mkdir(parents=True)
        config = Config(project_name=name, project_base=str(project_path))
        config.to_yaml(project_path)
        config.verify_config()

    @classmethod
    def create(cls, name: str, path: Path, config: Config) -> "Project":
        project_path = path / name
        project_path.mkdir(parents=True)
        config.verify_config()
        return cls(project_path, config) 

    @classmethod
    def load(cls, project_path: Path) -> "Project":
        config = Config.from_yaml(project_path / "config.yaml")
        config.verify_config()
        return cls(project_path, config)
    
    def make_database(self, threads=32, ani=0.98, c=200):
        """
        Makes the MAG-E database from the input genomes. 
        
        Each real sample is matched against the genomes of this database. Close
        genome matches are used to construct the mirror specification of the sample.
        """
        # make the directory housing the MAG-E database
        os.makedirs(self.config.maggie_db_dir)

        # the database table will be built up from the user-specified cluster assignments
        db_table = pd.read_csv(self.config.cluster_assignments, index=None)
        db.check_genome_files_exist(db_table.genome, self.config.genomes_dir)

        # Cluster the genomes at the strain level
        db.run_strain_clustering(db_table, self.config.genomes_dir, self.config.strain_dreps, threads, ani)

        # Add the strain cluster information to the database table
        db_table = db.build_database_table(db_table, self.config.genomes_dir, self.config.strain_dreps, self.config.maggie_db_md)

        # build syldb of species representatives
        reprs = db_table.FileLocation[db_table.isSpeciesRepr]
        db.construct_sylphdb(reprs, db_prefix='repr', t=threads,c=c,force=True)

        # build sylphdb of all genomes
        all_genomes = db_table.FileLocation
        db.construct_sylphdb(all_genomes, db_prefix='all', t=threads,c=c, force=True)

        return db_table
    
    def make_mirrors(self, threads, c):
        """
        Constructs the mirror specifications for each sample. 
        """

        # Make a Sylph sketch (paired.sylsp) of each sample.
        ms.construct_sylphsp(self.config.samples_dir, self.config.sylsp_dir, self.config.prefix1, t=threads,c=c)

        # Sylph query and profile each sample.
        ms.sylph_profile(self.config.sylsp_dir, self.config.maggie_db_dir, self.config.simulation_dir, t=threads)
        ms.sylph_query(self.config.sylsp_dir, self.config.maggie_db_dir, self.config.simulation_dir, t=threads)

        # Collect the Sylph results, and construct the mirror specifications.
        profiles = glob.glob(join(self.config.simulation_dir, '*_sylph_profiles.tsv'))
        queries = glob.glob(join(self.config.simulation_dir, '*_sylph_query.tsv'))
        for profile, query in zip(profiles, queries):
            ms.construct_metagenomic_specification(self.config.simulation_dir, self.config.maggie_db_md, profile, query)

    def simulate_mgx(self, threads, n_reads):
        """
        Simulated metagenomic data. 
        """
        specs = sorted(glob.glob(join(self.config.simulation_dir, '*_metagenomic.spec.csv')))
        for spec in specs:
            ms.run_InSilicoSeq(
                spec, self.config.simulation_dir, self.config.read_counts, 
                self.config.prefix1, n_reads, threads, force=True
            )
