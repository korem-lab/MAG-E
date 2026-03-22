import typer
from typer.core import TyperGroup
from pathlib import Path
from typing import Optional
from .project import Project

class OrderedGroup(TyperGroup):
    def list_commands(self, ctx):
        return list(self.commands)

app = typer.Typer(help='MAG-E -- Metagenome Assembled Genome Evaluator',cls=OrderedGroup)

@app.command()
def start_project(
    name: str = typer.Argument(..., help="New MAG-E project name."),
    directory: Path = typer.Argument(..., help="Project directory.")

):
    Project.create_empty(name, directory)
    "Creates new project with an empty configuration file."

@app.command()
def verify_config():
    "Makes sure everything needed is present in the configuration."
    pass

@app.command()
def get_current_project():
    "Prints the name and path of the project currently operating on."
    pass

@app.command()
def set_project(
    name: str = typer.Argument(..., help="MAG-E will operate on this project.")
):
    "Changes the current project to the one specified."
    pass

@app.command()
def make_maggie_db(
    directory: Path = typer.Argument(..., help="Project directory"),
    threads: int = typer.Option(32, help='Threads to run dRep clustering and sylph'),
    ani: float = typer.Option(0.98, help='ANI for strain clustering'),
    c: int = typer.Option(200, help='Sylph sketch density.')
):
    """
    Constructs the MAG-E database for the project.
    """
    project = Project.load(directory)
    project.make_database(threads, ani, c)

@app.command()
def make_mirrors(
    directory: Path = typer.Argument(..., help="Project directory"),
    threads: int = typer.Option(32, help='Threads to run dRep clustering and sylph'),
    c: int = typer.Option(200, help='Sylph sketch density.')
):
    """
    Constructs a mirror specification for each sample in the directory.
    """
    project = Project.load(directory)
    project.make_mirrors(threads, c)

@app.command()
def simulate_mgx(
    directory: Path = typer.Argument(..., help="Project directory"),
    threads: int = typer.Option(32, help='Threads to run InSilicoSeq'),
    n_reads: str = typer.Option('auto', help='Number of reads to simulate for each sample. Default matches the sample.')
):
    """
    Simulated metagenomes using the mirror specifications.
    """
    project = Project.load(directory)
    project.simulate_mgx(threads, n_reads)

@app.command()
def create_manifest():
    """
    Builds the manifest listing all MAG generation tasks to be run.
    """
    pass

@app.command()
def run_assembly():
    """
    Runs assembly tasks.
    """
    pass

@app.command()
def construct_ground_truth():
    """
    Constructs the ground truth for each assembly.
    """
    pass

@app.command()
def run_binning():
    """
    Runs binning and quality control tasks.
    """
    pass

@app.command()
def run_wrap():
    """
    Runs wrappers and quality control tasks.
    """
    pass

@app.command()
def calc_per_genome_metrics():
    """
    Calculate MAG-E recall, precision, F-score per genome. 
    """
    pass

@app.command()
def calc_contig_level_metrics():
    """
    Calculate MAG-E recall, precision, F-score of contig groups. 
    """
    pass

@app.command()
def construct_recoverable_set(
    precision = typer.Argument(..., help='Minimum MAG-E precision for the recoverable set.'),
    recall = typer.Argument(..., help='Minimum MAG-E recall for the recoverable set.')
):
    """
    Constructs the set of genomes considered recoverable. These will be used pipeline evaluations.
    """
    pass

@app.command()
def evaluate_pipelines(
    no_plots: bool = typer.Option(True, help='No plots will be produced.'),
    all_genomes: bool = typer.Option(False, help='Evaluate on all genomes rather than just the recoverable set.')
):
    """
    Builds a linear mixed model of the per-genome metrics over all MAG-pipelines.
    Estimated marginal means (and other values) from model are reported and plotted.
    """
    pass
