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
def initialize_project(
    name: str = typer.Argument(..., help="New MAG-E project name."),
    directory: Path = typer.Argument(..., help="Project directory.")

):
    Project.create_empty(name, directory)
    "Creates new project with an empty configuration file."

@app.command()
def build_project():
    "Verifies the config, and populates the project with the required directories."
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
    project.make_database(threads, ani, c, write_to_disk=True)


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
def construct_tasks(
    directory: Path = typer.Argument(..., help="Project directory"),
):
    """
    Constructs a task manifest, listing all MAG generation tasks to be run.
    """
    project = Project.load(directory)
    project.construct_tasks(write_to_disk=True)

@app.command()
def run_assembly(
    directory: Path = typer.Argument(..., help="Project directory"),
    threads: int = typer.Option(32, help='Threads to assemblers')
):
    """
    Runs assembly tasks.
    """
    project = Project.load(directory)
    project.run_assembly(threads)

@app.command()
def construct_ground_truth(
    directory: Path = typer.Argument(..., help="Project directory"),
    min_contig_len: int = typer.Option(100, help='Minimum length of contigs that can contribute to ground truth metrics.'), 
    min_pident: float = typer.Option(99, help='Minimum percent identity for a ground truth match'),
    min_aln_prop: float = typer.Option(99, help='alignment_length >= contig_length*(min_aln_prop) for a ground truth match'),
    max_aln_prop: float = typer.Option(101, help='alignment_length <= contig_length*(max_aln_prop) for a ground truth match')
):
    """
    Constructs the ground truth for each assembly.
    """
    project = Project.load(directory)
    project.construct_ground_truth(min_contig_len, min_pident, min_aln_prop, max_aln_prop)

@app.command()
def run_binning():
    """
    Runs binning and wrappers.
    """

@app.command()
def run_quality_control():
    """
    Runs quality control tools.
    """

@app.command()
def calc_per_genome_metrics():
    """
    Calculate MAG-E recall, precision, F-score per genome. 
    """

@app.command()
def calc_contig_level_metrics():
    """
    Calculate MAG-E recall, precision, F-score of contig groups. 
    """

@app.command()
def construct_recoverable_set(
    precision = typer.Argument(..., help='Minimum MAG-E precision for the recoverable set.'),
    recall = typer.Argument(..., help='Minimum MAG-E recall for the recoverable set.')
):
    """
    Constructs the set of genomes considered recoverable. These will be used pipeline evaluations.
    """

@app.command()
def evaluate_pipelines(
    no_plots: bool = typer.Option(True, help='No plots will be produced.'),
    all_genomes: bool = typer.Option(False, help='Evaluate on all genomes rather than just the recoverable set.')
):
    """
    Builds a linear mixed model of the per-genome metrics over all MAG-pipelines.
    Estimated marginal means (and other values) from model are reported and plotted.
    """
