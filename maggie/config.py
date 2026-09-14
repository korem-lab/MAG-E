import os
import glob
from os.path import join, exists
from dataclasses import dataclass, field, asdict, fields
from collections import defaultdict
from pathlib import Path
import yaml
import pandas as pd
from .utils import tupleize
from .tasks import quality_control, assembly, binning


def represent_list(dumper, data):
    return dumper.represent_sequence('tag:yaml.org,2002:seq', data, flow_style=True)
yaml.add_representer(list, represent_list)

def make_dict(data):
    data=tupleize(data)
    ret=defaultdict(set)
    for t,p in data:
        ret[t].add(p)
    return ret

null = 'NULL'
null_asm = [['ASM', 'OPT'], ['ASM', 'OPT']]
null_bin = [['BIN', 'OPT'], ['BIN', 'OPT']]
null_cov = [['COV', 'OPT'], ['COV', 'OPT']]
null_mode = ['NULL', 'NULL']
null_qc = ['NULL', 'NULL']
null_refiner = [
    ['REFN', 'REFN_OPT', 'ASM','ASMOPT', 'COV', 'COVOPT', [['BIN', 'BINOPT', 'MODE'], ['BIN', 'BINOPT', 'MODE']]]
]

@dataclass
class Config:
    project_name: str = null
    project_base: str = null
    use_api: str = null
    samples_dir: str = null
    read_counts: str = null
    contig_properties: str = null
    prefix1: str  = null
    ecosystem_db: str = null
    assemblers: list = field(default_factory=lambda: null_asm)
    binners: list = field(default_factory=lambda: null_bin)
    coverage: list = field(default_factory=lambda: null_cov)
    binning_modes: list =  field(default_factory=lambda: null_mode)
    refiners: list = field(default_factory=lambda: null_refiner)
    qctools: list = field(default_factory=lambda: null_qc)
    simulate: str = 'yes'

    @classmethod
    def from_yaml(cls, path: Path) -> "Config":
        with open(path) as f:
            data = yaml.safe_load(f)
        keys = {f.name for f in fields(cls)}
        for k in keys:
            if k not in data:
                data[k] = []
        o = cls(**data)
        o.read_counts = join(o.project_base, o.read_counts)
        return o

    def to_yaml(self, path: Path) -> None:
        with open(join(path, 'config.yaml'), "w") as f:
            yaml.dump(asdict(self), f, default_flow_style=False)

    def verify_config(self):
        """Check whether the configuration file is valid. """
        assert exists(self.project_base)
        assert exists(self.samples_dir)
        assert exists(self.read_counts)
        # Either left blank, or set to no or yes
        assert (self.simulate == [] or self.simulate == 'no' or self.simulate == 'yes')
        # If we're simulating, we need to have a valid ecosystem database
        if (self.simulate == 'yes'):
            assert exists(self.ecosystem_db)
            assert exists(self.ecosystem_db_metadata)
            assert exists(self.genomes_dir)
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
        self.coverage = [tupleize(e) for e in self.coverage]

    def get_list_of(self, k):
        ret = ''
        if k == 'targets':
            targets = [os.path.basename(e) for e in glob.glob(f'{self.simulation_dir}/*.fastq.gz')]
            targets = [e.split(f'_{self.prefix1}')[0] for e in targets if f'_{self.prefix1}' in e]
            ret = ' '.join(targets)
        elif k == 'assemblers':
            ret = ' '.join(set(a for a,p in self.assemblers))
        elif k == 'binners':
            ret = ' '.join(set(b for b,p in self.binners))
        elif k == 'coverage':
            ret = ' '.join(set(b for b,p in self.coverage))
        elif k == 'modes':
            ret = ' '.join(set(self.binning_modes))
        elif k == 'qctools':
            ret = ' '.join(set(self.qctools))
        elif k == 'refiners':
            ret = ' '.join(set(r for r,*rest in self.refiners))
        else:
            Exception('--of is invalid.')
        return ret

    def get_options_within(self, tool):
        ret = ''
        if tool in [a for a, p in self.assemblers]:
            ret = '\0'.join(p for a, p in self.assemblers if a == tool)
        elif tool in [a for a, p in self.binners]:
            ret = '\0'.join(p for a, p in self.binners if a == tool)
        elif tool in [a for a, p in self.coverage]:
            ret = '\0'.join(p for a, p in self.coverage if a == tool)
        elif tool in [a for a, *rest in self.refiners if a == tool]:
            ret = '\0'.join(p for a,p, *rest in self.refiners if a == tool)
        else:
            Exception(f'{tool} is invalid.')
        return ret

    def get_binner_sets_within(
        self, refiner, refiner_option, assembler, 
        assembler_option, coverage, coverage_option, **kwargs
    ):
        binner_sets = list()
        for rf, ro, a, ao, c, co, binner_set in self.refiners:
            if refiner == rf and ro == refiner_option and assembler == a and assembler_option == ao and coverage == c and coverage_option == co:
                binner_sets.append(':'.join([','.join(b) for b in binner_set]))
        return '\0'.join(binner_sets)

    def get_samples_within_mode(self, mode, target):
        assert mode in self.binning_modes, Exception("Supplied binning mode is invalid.")
        df = pd.read_csv(os.path.join(self.project_base, f'{mode}_datasets.csv'))
        return ' '.join(df.loc[df.target == target].dataset.values)

    @property
    def ecosystem_db_metadata(self):
        return join(self.ecosystem_db, 'ecosystem_metadata.csv')

    @property
    def genomes_dir(self):
        return join(self.ecosystem_db, 'genomes')

    @property
    def strain_dreps(self):
        return join(self.ecosystem_db, 'strain_dreps')

    @property
    def simulation_dir(self):
        if self.simulate == [] or self.simulate == 'yes':
            return join(self.project_base, 'simulations')
        else:
            return join(self.project_base, 'read_links')

    @property
    def sylsp_dir(self):
        return join(self.simulation_dir, 'sylsp')

    @property
    def assembly_cache(self):
        return join(self.project_base, 'assembly_cache')

    @property
    def coverage_cache(self):
        return join(self.project_base, 'coverage_cache')

    @property
    def bintask_dir(self):
        return join(self.project_base, 'bintask_dir')

    @property
    def manifest(self):
        return join(self.project_base, 'manifest.csv')

    @property
    def evaluation_dir(self):
        return join(self.project_base, 'evaluation')