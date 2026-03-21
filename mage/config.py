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
    prefix: str  = 'NULL'
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
        all_good &= os.path.exists(Config.project_name)
        all_good &= os.path.exists(Config.genomes_dir)
        all_good &= os.path.exists(Config.cluster_assignments)
        all_good &= os.path.exists(Config.samples_dir)
        all_good &= os.path.exists(Config.read_counts)
        all_good &= os.path.exists(Config.contig_properties)
        all_good &= 'NULL' not in Config.prefix
        all_good &= Config.n_samples > 0
        print('Config looks good.' if all_good else 'Config misspecified. Consider revising.')

    @property
    def mage_db_dir():
        return join(Config.project_base, 'database')

    @property
    def mage_db_md():
        return join(Config.project_base, 'database/metadata.csv')

    @property
    def strain_dreps():
        return join(Config.project_base, 'database/strain_dreps')
