import os
from os.path import join, exists
from dataclasses import dataclass, field, asdict
from pathlib import Path
import yaml
from .utils import tupleize
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
    read_counts: str = 'NULL'
    contig_properties: str = 'NULL'
    prefix1: str  = 'NULL'
    assemblers: list = field(default_factory=lambda: [['ASM', 'OPT'], ['ASM', 'OPT']])
    binners: list = field(default_factory=lambda: [['BIN', 'OPT'], ['BIN', 'OPT']])
    binning_modes: list =  field(default_factory=lambda: ['MODE', 'MODE'])
    refiners: list = field(default_factory=lambda: 
        [
           ['REFN', 'REFN_OPT', [['ASM', 'ASMOPT', 'BIN', 'BINOPT', 'MODE'], ['ASM', 'ASMOPT', 'BIN', 'BINOPT', 'MODE']]]
        ]
    )
    qctools: list = field(default_factory=lambda: ['NULL', 'NULL'])

    @classmethod
    def from_yaml(cls, path: Path) -> "Config":
        with open(path) as f:
            data = yaml.safe_load(f)
        o = cls(**data)
        o.cluster_assignments = join(o.project_base, o.cluster_assignments)
        o.read_counts = join(o.project_base, o.read_counts)
        return o

    def to_yaml(self, path: Path) -> None:
        with open(join(path, 'config.yaml'), "w") as f:
            yaml.dump(asdict(self), f, default_flow_style=False)

    def verify_config(self):
        """Check whether the configuration file is valid. """
        assert exists(self.project_base)
        assert exists(self.genomes_dir)
        assert exists(self.samples_dir)
        assert exists(self.cluster_assignments)
        assert exists(self.read_counts)
        assert 'NULL' not in self.prefix1
        assert all(exists(join(self.project_base, f'{mode}_datasets.csv')) for mode in self.binning_modes)

        # If the api is being used, then make sure every assembly, bin, and refiner class is implemented
        if self.use_api == 'yes':
            for ab,_ in self.assemblers:
                assert hasattr(assembly, ab+'Assembler')
            for bn,_ in self.binners:
                assert(hasattr(binning, bn+'Binner'))
            for rf, _, pipeline in self.refiners:
                assert(hasattr(binning, rf+'Refiner'))
                for (ab, _, bn, _, _) in pipeline:
                    assert(hasattr(assembly, ab+'Assembler'))
                    assert(hasattr(binning, bn+'Binner'))
        
        # Make sure the quality control API is implemented 
        for q in self.qctools:
            assert hasattr(quality_control, q+'QCTool')
        
        self.assemblers = [tupleize(e) for e in self.assemblers]
        self.binners = [tupleize(e) for e in self.binners]
        self.refiners = [tupleize(e) for e in self.refiners]

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
    def simulation_dir(self):
        return join(self.project_base, 'simulations')

    @property
    def sylsp_dir(self):
        return join(self.simulation_dir, 'sylsp')

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