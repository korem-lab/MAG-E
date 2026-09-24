import pandas as pd
from os.path import exists, join
from itertools import product
import hashlib
import ast
import json
from .utils import flatten

cols = [
    'assembler', 'assembler_options', 'coverage', 'coverage_options', 
    'binner', 'binner_options', 'binning_mode',
    'assembler_dir', 'coverage_dir', 'binner_dir', 'simulation_dir', 'qctools', 
    'assembler_summary_name', 'binner_summary_name', 
    'coverage_summary_name', 'target', 'samples', 'task_hash', 
    'assembler_task_id', 'coverage_task_id', 'binner_task_id'
]

def inc(n):
    c=n+1 if pd.notna(n) else 1
    while True:
        yield c
        c+=1

def manifest_format(rf, ro, a, ao, c, co, binner_set, target):
    b = [e[0] for e in binner_set]
    bo = [e[1] for e in binner_set]
    mode = [e[2] for e in binner_set]
    return (a, ao, c, co, [rf]+[b], [ro]+[bo], [None]+[mode], target)

def fill(df, taskcols, idcol):
    # get the unique set of tasks
    tasks = df[taskcols + [idcol]].drop_duplicates()
    tasks = tasks.sort_values(idcol).drop_duplicates(subset=taskcols, keep='first')
    # assign every unique task a unique id, preserving already assigned ids
    _id = inc(tasks[idcol].max())
    tasks[idcol] = tasks[idcol].apply(lambda x: x if not pd.isna(x) else next(_id))
    # make this a lookup on the taskcols
    tasks.set_index(taskcols,inplace=True)
    # remap the ids
    return pd.Series(list(map(tasks[idcol].get, zip(*[df[e] for e in taskcols]))), index=df.index).astype(int)


class Manifest():
    def __init__(self, project_base):
        self.manifest_fl = join(project_base, '_manifest.csv')
        if exists(self.manifest_fl):
            self.read()
        else:
            self.m = pd.DataFrame(columns=cols)
        self.check()

    def check(self):
        assert(self.m.notna().all().all())
        assert(not self.m.empty)

    def sync(self, cfg):
        prodcols = [
            'assembler', 'assembler_options', 'coverage', 'coverage_options', 
            'binner', 'binner_options', 'binning_mode', 'target'
        ]
        targets = sorted(cfg.get_list_of('targets').split())
        _m = pd.DataFrame(
            [flatten(e) for e in product(cfg.assemblers, cfg.coverage, cfg.binners, cfg.binning_modes, targets)], 
            columns=prodcols
        )
        _rm = pd.DataFrame([
            manifest_format(r,t) for r in cfg.refiners for t in targets], 
            columns=prodcols
        )
        _m = pd.concat([_m, _rm])
        self.m = pd.merge(self.m, _m, how='outer', on=prodcols)

        # fill task ids
        self.fill(['assembler', 'assembler_options'], 'assembler_task_id')
        self.fill(['assembler', 'assembler_options', 'coverage', 'coverage_options'], 'coverage_task_id')
        self.fill([
            'assembler', 'assembler_options', 'coverage', 
            'coverage_options', 'binner', 'binner_options', 'binning_mode'], 'binner_task_id'
        )
        self.m = self.m.sort_values(
            ['assembler_task_id', 'coverage_task_id', 'binner_task_id']
        )

        # fill summary names
        self.fill_sname('assembler')
        self.fill_sname('coverage')
        self.fill_sname('binner')

        # constants
        self.m['simulation_dir'] = cfg.simulation_dir
        self.m['qctools'] = [tuple(cfg.qctools)] * len(self.m)

        # Constructs the paths where each particular task output (e.g MEGAHIT assembly task)
        # will be written within the project
        self.set_task_locations('assembler', cfg.project_base)
        self.set_task_locations('coverage', cfg.project_base)
        self.set_task_locations('binner', cfg.project_base)

        self.add_sample_list(cfg)
        self.add_task_hash()

        self.check()
        self.write()

    def read(self):
        df = pd.read_csv(self.manifest_fl, index_col=0)
        df['samples'] = df.samples.apply(ast.literal_eval)
        if 'qctools' in df.columns:
            df.qctools = df.qctools.apply(ast.literal_eval)
        self.m = df

    def write(self):
        self.m.to_csv(self.manifest_fl)

    def add_task_hash(self):
        cols = list(set(self.m.columns) - {'task_hash'})
        self.m['task_hash'] = self.m.apply(lambda x: hashlib.sha256(json.dumps(x[cols].to_dict(), sort_keys=True).encode()).hexdigest(), axis=1)

    def fill(self, taskcols, idcol):
        self.m[idcol] = fill(self.m, taskcols, idcol)

    def fill_sname(self, tool):
        opts = f'{tool}_options'
        sname = f'{tool}_summary_name'
        parts = [fill(g, [tool, opts], sname) for _,g in self.m.groupby(tool)]
        res = pd.concat(parts)
        self.m[sname] = res.reindex(self.m.index)

    def set_task_locations(self, tool, basedir):
        self.m[f'{tool}_dir'] = self.m.apply(
            lambda x: join(basedir, f'{tool}_tasks', 'task_'+str(x[f'{tool}_task_id']), str(x.target)), axis=1
        )

    def add_sample_list(self, cfg):
        def get_samples(m,t):
            samples = cfg.get_samples_within_mode(m,t).split()
            return [t] + sorted(list(set(samples) - {t}))
        pairs = self.m[['binning_mode', 'target']].drop_duplicates()
        pairs['samples'] = pairs.apply(lambda x: get_samples(x.binning_mode, x.target),axis=1)
        pairs.set_index(['binning_mode', 'target'], inplace=True)
        self.m['samples'] = list(map(pairs['samples'].get, zip(self.m['binning_mode'], self.m['target'])))

    def get_tasks(self, 
            assembler, aopt, coverage, copt, binner, bopt, binning_mode, binner_set, target, cov_sample
        ):
        flt = lambda x,y: self.m[x] == y if y is not None else pd.Series(True, index=self.m.index)
        ms_flt = lambda ms: self.m.samples.apply(lambda x: ms in x) if ms is not None else pd.Series(True, index=manifest.index)
        if binner_set:
            binner_set = tuple(tuple(e.split(',')) for e in binner_set.split(':'))
        return self.m.loc[
            flt('target', target) & flt('assembler', assembler) & flt('assembler_options', aopt) &
            flt('coverage', coverage) & flt('coverage_options', copt) & flt('binner', binner) & 
            flt('binner_options', bopt) & flt('binning_mode', binning_mode) &
            flt('binner_set', binner_set) & ms_flt(cov_sample)
        ]

    def get_task_directory(self, **kwargs):
        tsks = self.get_tasks(cov_sample=None, **kwargs)
        assembler, coverage, binner = kwargs['assembler'], kwargs['coverage'], kwargs['binner']
        binning_mode, target = kwargs['binning_mode'], kwargs['target']
        if assembler and coverage and binner and binning_mode and target:
            tsks = tsks[['binner_dir']].drop_duplicates()
        elif assembler and coverage and target:
            tsks = tsks[['coverage_dir']].drop_duplicates()
        elif assembler and target:
            tsks = tsks[['assembly_dir']].drop_duplicates()
        else:
            Exception("Invalid query.")
        if tsks.empty:
            raise Exception("No tasks fit the query.")
        elif len(tsks) > 1:
            raise Exception("More than task fits the query") 
        return tsks.iloc[0,0]

