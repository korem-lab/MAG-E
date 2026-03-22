import os
from os.path import join
from dataclasses import dataclass, field, asdict
from pathlib import Path
import yaml

@dataclass
class Config:
    project_name: str = 'NULL'
    project_base: str = 'NULL'
    genomes_dir: str = 'NULL'
    cluster_assignments: str = 'NULL'
    samples_dir: str = 'NULL'
    read_counts: str = 'NULL'
    contig_properties: str = 'NULL'
    prefix1: str  = 'NULL'
    n_samples: int = -1

    @classmethod
    def from_yaml(cls, path: Path) -> "Config":
        with open(path) as f:
            data = yaml.safe_load(f)
        return cls(**data)

    def to_yaml(self, path: Path) -> None:
        with open(path / 'config.yaml', "w") as f:
            yaml.dump(asdict(self), f, default_flow_style=False)

    def verify_config(self):
        """Check various properties about the project to make sure its all good. """
        all_good = True
        all_good &= os.path.exists(self.project_name)
        all_good &= os.path.exists(self.genomes_dir)
        all_good &= os.path.exists(self.cluster_assignments)
        all_good &= os.path.exists(self.samples_dir)
        all_good &= os.path.exists(self.read_counts)
        all_good &= os.path.exists(self.contig_properties)
        all_good &= 'NULL' not in self.prefix
        all_good &= self.n_samples > 0
        print('Config looks good.' if all_good else 'Config misspecified. Consider revising.')

    @property
    def maggie_db_dir(self):
        return join(self.project_base, 'database')

    @property
    def maggie_db_md(self):
        return join(self.maggie_db_dir, 'metadata.csv')

    @property
    def strain_dreps(self):
        return join(self.maggie_db_dir, 'strain_dreps')

    @property
    def sylsp_dir(self):
        return join(self.project_base, 'sylsp')

    @property
    def simulation_dir(self):
        return join(self.project_base, 'simulations')
