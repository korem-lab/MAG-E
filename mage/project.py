import pandas as pd
import os
from os.path import join, exists
from pathlib import Path
from .config import Config

import database as db


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
    
    def make_database(self, threads=32, ani=0.98):
        """
        Makes the MAG-E database from the input genomes. 
        
        Each real sample is matched against the genomes of this database. Close
        genome matches are used to construct the mirror specification of the sample.
        """
        # make the directory housing the MAG-E database
        os.makedirs(self.config.mage_db_dir)

        db_table = pd.read_csv(self.config.cluster_assignments, index=None)
        db.check_each_genome_file_exists(db_table.genome, self.config.genomes_dir)

        # Cluster the genomes at the strain level
        db.run_strain_dreps(db_table, self.config.genomes_dir, self.config.strain_dreps, threads, ani)

        # metadata for database
        db.build_database_table(db_table, self.config.genomes_dir, self.config.strain_dreps, self.config.mage_db_md)

        # build syldb of species representatives
        reprs= db_table.FileLocation[db_table.isSpeciesRepr]
        db.construct_sylphdb(reprs, db_prefix='repr', t=threads,force=True)

        # build sylphdb of all genomes
        all_genomes = db_table.FileLocation
        db.construct_sylphdb(all_genomes, db_prefix='all', t=threads,force=True)
