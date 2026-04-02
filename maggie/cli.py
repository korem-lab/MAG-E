import typer
from os.path import exists, join
from os import remove
import pandas as pd
from typer.core import TyperGroup
from pathlib import Path
from typing import Optional
from .project import Project

class OrderedGroup(TyperGroup):
    def list_commands(self, ctx):
        return list(self.commands)

app = typer.Typer(help='MAG-E -- Metagenome Assembled Genome Evaluator',cls=OrderedGroup)

def add_new_project_to_cache(name, directory):
    MAGE_dir = Path(__file__).parent.parent
    cache_name = join(MAGE_dir, '.project_cache.tsv')
    record = {'name':name, 'directory': directory, 'is_current':False}
    if not exists(cache_name):
        # if its the first project in existence, set is current to true
        record['is_current'] = True
        df = pd.DataFrame([record])
    else:
        df = pd.read_csv(cache_name, sep='\t')
        assert df[(df.name == name)].empty, Exception(f'Project already exists with name {name}, pick another name.')
        df.loc[len(df)] = [name, directory, False]
    df.to_csv(cache_name, sep='\t', index=None)

def set_project_to_current(name):
    MAGE_dir = Path(__file__).parent.parent
    cache_name = join(MAGE_dir, '.project_cache.tsv')
    df = pd.read_csv(cache_name, sep='\t')
    assert df.is_current.sum() == 1 # exactly one current allowed
    df.is_current = False # reset

    # Set the new current
    assert name in df.name.to_list(), Exception(f'{name} not in cached projects.')
    current_idx = df.index[df.name == name]
    df.loc[current_idx.item(), 'is_current'] = True
    df.to_csv(cache_name, sep='\t', index=None)

def get_current_project():
    MAGE_dir = Path(__file__).parent.parent
    cache_name = join(MAGE_dir, '.project_cache.tsv')
    df = pd.read_csv(cache_name, sep='\t')
    return df[df.is_current].name.item(), df[df.is_current].directory.item()

@app.command()
def initialize_project(
    directory: Path = typer.Argument(..., help="Project directory."),
    name: str = typer.Argument(..., help="New MAG-E project name."),
    use_api: str = typer.Option('yes', help="Will use API classes only (yes), not at all (no), or when available (when_possible)")
):
    "Creates new project with an empty configuration file."
    assert use_api in ['yes', 'no', 'when_possible']
    Project.create_empty(name, directory, use_api)
    add_new_project_to_cache(name, directory/name)
    set_project_to_current(name)

@app.command()
def get_current():
    "Prints the name and path of the project currently operating on."
    name, dir = get_current_project()
    print(f'{name}\t{dir}', flush=True)

@app.command()
def set_current(
    name: str = typer.Argument(..., help="MAG-E will operate on this project.")
):
    "Changes the current project to the one specified."
    set_project_to_current(name)

@app.command()
def reset_cache():
    remove(join(Path(__file__).parent.parent,'.project_cache.tsv'))

@app.command()
def verify_project(
):
    "Verifies the config, and populates the project with the required directories."
    _, directory = get_current_project()
    project = Project.load(directory, verify=False)
    project.build_core_directories()
    project.config.verify_config()

@app.command()
def make_maggie_db(
    threads: int = typer.Option(8, help='Threads to run dRep clustering and sylph'),
    ani: float = typer.Option(0.98, help='ANI for strain clustering'),
    c: int = typer.Option(200, help='Sylph sketch density.')
):
    """
    Constructs the MAG-E database for the project.
    """
    _, directory = get_current_project()
    project = Project.load(directory)
    project.make_database(threads, ani, c, write_to_disk=True)


@app.command()
def make_mirrors(
    threads: int = typer.Option(32, help='Threads to run dRep clustering and sylph'),
    c: int = typer.Option(200, help='Sylph sketch density.')
):
    """
    Constructs a mirror specification for each sample in the directory.
    """
    _, directory = get_current_project()
    project = Project.load(directory)
    project.make_mirrors(threads, c)

@app.command()
def simulate_mgx(
    threads: int = typer.Option(32, help='Threads to run InSilicoSeq'),
    n_reads: str = typer.Option('auto', help='Number of reads to simulate for each sample. Default matches the sample.')
):
    """
    Simulated metagenomes using the mirror specifications.
    """
    _, directory = get_current_project()
    project = Project.load(directory)
    project.simulate_mgx(threads, n_reads)

@app.command()
def construct_tasks(
):
    """
    Constructs a task manifest, listing all MAG generation tasks to be run.
    """
    _, directory = get_current_project()
    project = Project.load(directory)
    project.construct_tasks(write_to_disk=True)

@app.command()
def run_assembly(
    threads: int = typer.Option(32, help='Threads to assemblers')
):
    """
    Runs assembly tasks.
    """
    _, directory = get_current_project()
    project = Project.load(directory)
    project.run_assembly(threads)

@app.command()
def construct_ground_truth(
    min_contig_len: int = typer.Option(100, help='Minimum length of contigs that can contribute to ground truth metrics.'), 
    min_pident: float = typer.Option(99, help='Minimum percent identity for a ground truth match'),
    min_aln_prop: float = typer.Option(99, help='alignment_length >= contig_length*(min_aln_prop) for a ground truth match'),
    max_aln_prop: float = typer.Option(101, help='alignment_length <= contig_length*(max_aln_prop) for a ground truth match')
):
    """
    Constructs the ground truth for each assembly.
    """
    _, directory = get_current_project()
    project = Project.load(directory)
    project.construct_ground_truth(min_contig_len, min_pident, min_aln_prop, max_aln_prop)

@app.command()
def run_binning(
    threads: int = typer.Argument(..., help='Number of threads to run each binning task.'),
    run_only: str = typer.Option(None, help='Run prep or binning only.'),
    force_bin_prep: bool = typer.Option(False, help='Will force rerun the prep.'),
    force_bin: bool = typer.Option(False, help='Will force rerun the binning.'),
    force_refine: bool = typer.Option(False, help='Will force rerun refining.')
):
    """
    Runs binning and refiners.
    """
    _, directory = get_current_project()
    project = Project.load(directory)
    project.run_binning(run_only, force_bin_prep, force_bin, threads)
    project.run_refine(run_only, force_refine, threads)

@app.command()
def run_quality_control(
    threads: int = typer.Argument(..., help='Number of threads to run each binning task.'),
):
    """
    Runs quality control tools.
    """
    _, directory = get_current_project()
    project = Project.load(directory)
    project.run_quality_control(threads)

@app.command()
def calc_per_genome_metrics(
):
    """
    Calculate MAG-E recall, precision, F-score per genome. 
    """
    _, directory = get_current_project()
    project = Project.load(directory)
    project.calc_per_genome_metrics()

@app.command()
def calc_contig_level_metrics():
    """
    Calculate MAG-E recall, precision, F-score of contig groups. 
    """

@app.command()
def evaluate_pipelines(
    no_plots: bool = typer.Option(True, help='No plots will be produced.'),
    precision: float = typer.Option(0.9, help='Minimum MAG-E precision for the recoverable set.'),
    recall: float = typer.Option(0.7, help='Minimum MAG-E recall for the recoverable set.')
):
    """
    Builds a linear mixed model of the per-genome metrics over all MAG-pipelines.
    Estimated marginal means (and other values) from model are reported and plotted.
    """
    _, directory = get_current_project()
    project = Project.load(directory)
    project.evaluate_pipelines(precision, recall, no_plots)


@app.command()
def compare_quality_control(
    no_plots: bool = typer.Option(True, help='No plots will be produced.'),
    precision = typer.Argument(..., help='Minimum MAG-E precision for the recoverable set.'),
    recall = typer.Argument(..., help='Minimum MAG-E recall for the recoverable set.')
):
    """
    Builds a linear mixed model of the per-genome metrics over all MAG-pipelines.
    Estimated marginal means (and other values) from model are reported and plotted.
    """