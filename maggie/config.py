import os
from os.path import join, exists
from dataclasses import dataclass, field, asdict
from pathlib import Path
import yaml
from .tasks import quality_control, assembly, binning

def represent_list(dumper, data):
    return dumper.represent_sequence('tag:yaml.org,2002:seq', data, flow_style=True)
yaml.add_representer(list, represent_list)

@dataclass
class Config:
    project_name: str = 'NULL'
    project_base: str = 'NULL'
    use_api: str = 'NULL'
    genomes_dir: str = 'NULL'
    cluster_assignments: str = 'NULL'
    samples_dir: str = 'NULL'
    read_counts: str = 'NULL' # you have to place this in the project directory after its constructed 'project_base/read_counts.csv'
    contig_properties: str = 'NULL'
    prefix1: str  = 'NULL'
    assemblers: list = field(default_factory=lambda: [('ASM', 'OPT'), ('ASM', 'OPT')])
    binners: list = field(default_factory=lambda: [('BIN', 'OPT'), ('BIN', 'OPT')])
    binning_modes: list =  field(default_factory=lambda: ['MODE', 'MODE'])# for example ['single', 'all', 'mash20', 'mash5']. Need files 'project_base/all_datasets.csv', 'project_base/mash20_datasets.csv' that you need to put there
    refiners: list = field(default_factory=lambda: [('REFN', (('ASM', 'ASMOPT', 'BIN', 'BINOPT', 'MODE'), ('ASM', 'ASMOPT', 'BIN', 'BINOPT', 'MODE')))])
    qctools: list = field(default_factory=lambda: ['NULL', 'NULL'])

    @classmethod
    def from_yaml(cls, path: Path) -> "Config":
        with open(path) as f:
            data = yaml.safe_load(f)
        return cls(**data)

    def to_yaml(self, path: Path) -> None:
        with open(path / 'config.yaml', "w") as f:
            yaml.dump(asdict(self), f, default_flow_style=False)

    def verify_config(self):
        """Check whether the configuration file is valid. """
        all_good = True
        all_good &= exists(self.project_base)
        all_good &= exists(self.genomes_dir)
        all_good &= exists(self.cluster_assignments)
        all_good &= exists(self.samples_dir)
        all_good &= exists(join(self.project_base, self.read_counts))
        all_good &= 'NULL' not in self.prefix1
        all_good &= all(exists(join(self.project_base, f'{mode}_datasets.csv') for mode in self.binning_modes))

        # If the api is being used, then make sure every assembly, bin, and refiner class is implemented
        if self.use_api == 'yes':
            for ab,_ in self.assemblers:
                assert hasattr(assembly, ab+'Assembler')
            for bn,_ in self.binners:
                assert(hasattr(binning, bn+'Binner'))
            for rf, pipeline in self.refiners:
                assert(hasattr(binning, rf+'Refiner'))
                for (ab, _, bn, _, _) in pipeline:
                    assert(hasattr(assembly, ab+'Assembler'))
                    assert(hasattr(binning, bn+'Assembler'))
        
        # Make sure the quality control API is implemented 
        for q in self.qctools:
            assert hasattr(quality_control, q+'QCTool')

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

    @property
    def assembly_cache(self):
        return join(self.project_base, 'assembly_cache')

    @property
    def bintask_dir(self):
        return join(self.project_base, 'bintask_dir')

    @property
    def manifest(self):
        return join(self.project_base, 'manifest.csv')

    @property
    def evaluation_dir(self):
        return join(self.project_base, 'evaluation')