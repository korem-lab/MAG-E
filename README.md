# MAG-E -- Metagenome Assembled Genome Evaluator

Welcome to MAG-E, a framework for evaluating metagenome assembled genome (MAG) construction pipelines at scale with ground truth!

This readme provides an explaination of how to run MAG-E on real metagenomic datasets to: 
- Construct realistic simulations of the dataset
- Evaluate MAG-generation pipelines on the simulated dataset dataset
## Installation
We gotta pip install InSilicoSeq. We need pandas==2.3.3 to work with dRep. 
## Setting up new projects

MAG-E operates over project directories, which have a specific structure and contents.
We begin by initializing a new project.

### Initializing a project
```
python -m maggie initalize-project  /absolute/path/to/maggie/projects/ myproject
```

Running `initialize-project` constructs a new project `myproject` at `/absolute/path/to/maggie/projects`
and, inside the project directory, writes a configuration yaml file: 
```
(maggie) izaak@dhcp-10-118-17-29 MAG-E % ls ../testingmaggie/myproject 
config.yaml
```

### Configuring MAG-E

MAG-E uses the `config.yaml` to organize the entire project. 
The inputs MAG-E requires and the pipelines MAG-E will evaluate are both specified through the `config.yaml`. After initialization, `config.yaml` looks like: 

```
project_base: /Users/izaak/Library/CloudStorage/Dropbox/Documents/Projects/testingmaggie/myproject
project_name: myproject
use_api: 'yes'
genomes_dir: 'NULL'
samples_dir: 'NULL'
prefix1: 'NULL'
read_counts: 'NULL'
cluster_assignments: 'NULL'
assemblers: [[ASM, OPT], [ASM, OPT]]
binners: [[BIN, OPT], [BIN, OPT]]
binning_modes: [MODE, MODE]
refiners: [[REFN, [[ASM, ASMOPT, BIN, BINOPT, MODE], [ASM, ASMOPT, BIN, BINOPT, MODE]]]]
qctools: ['NULL', 'NULL']
```

Three fields are already populated: the project base path, name, and the `use_api` token set to `yes`, which determines how the bins from each MAG-generation pipeline will be provided to MAG-E (explained further below).

We must now populate the inputs required by MAG-E:
- `genomes_dir`:
Absolute path to directory containing all genomes that will form the MAG-E database. Genomes must be gzipped fasta files (prefix `.fasta.gz`).

- `samples_dir`:
Absolute path to directory containing all gzipped fastq files (prefix `.fastq.gz`).

- `prefix1`:
prefix of pair 1, typically either `_R1` or `_1`

- `read_counts`:
Name of the csv file specifying the number of read pairs for each sample. This must be moved
inside the project directory.
It should have the format:
```
sample,count
SampleA,10000
...
```

- `cluster_assignments`:
The name of the csv, specifying the species-level cluster assignments of each genome in `genomes_dir`. 
This csv must be moved into the project directory.
It has four required fields: `genome`, `SpeciesRepr`, and `GenomeType` and `N50`, `Length`. `genome` is the name 
of a genome which must exist at `genomes_dir/<genome>.fasta.gz`. `SpeciesRepr`
specifies the species-level cluster assignments. Genomes with the `SpeciesRepr` value are in the same cluster, the `SpeciesRepr` value (which must be a genome `genomes`) acts as the cluster representative, `N50` is the N50 of each genome fasta.
`GenomeType` specifies whether the genome is a metagenome-assembled genome, or sequenced as an isolate. 
Note that, MAG-E will only evaluate performance against the Isolate genomes. An example of the structure for the `cluster_assignments` csv is as follows:
```
genome,SpeciesRepr,GenomeType
genomeA,genomeA,Isolate
genomeB,genomeA,MAG
genomeC,genomeA,Isolate
genomeD,genomeD,Isolate
```

We must also specify the MAG-generation pipelines MAG-E will evaluate.
A MAG-generation pipeline is a combination of an assembler (and variable options), binner (and variable options), and binning mode. Pipelines are specified by populating the following:
- `assemblers`:
`assemblers` must be a list of lists. Each element of the outer-list is a 2-tuple, consisting of the name
of an assembler and the command line options used to run the assembler. `default` is used
to run the assembler with default options. For example, say we want to evaluate pipelines with MEGAHIT
in default, MEGAHIT with no bubble merging but --min-count of 5, and metaSPAdes in default. We would populate `assemblers` with `[['MEGAHIT', 'default'], ['MEGAHIT', '--no-bubble --min-count 5'], ['metaSPAdes', 'default]]`.

- `binners`:
`binners` follows the same format as `assemblers`. So, `[['METABAT2','default'], ['CONCOCT','default']]` would evaluate METABAT2 and CONCOCT MAG-pipelines in default. 

- `binning_mode`:
Binning mode determines which sample's information will be used to help bin a target sample. For example,
in `single` only information from the target sample will be used. Binning mode is very general, 
and allows any set of user-defined samples to be used to bin bin the target. For example, say,
we want to try `single`, `all`, `mash20`, and `mash3`, where `all` uses all samples in the dataset.
`mash20` uses the 20 closest samples to the target defined by mash distance, 
and `mash3`, which uses the closest 3. To achieve this, you must populate
binning mode with `['single', 'all', 'mash20', 'mash3']` and in the project directory, place
the files `all_datasets.csv`, `mash20_datasets.csv`, `mash3_datasets.csv`. 
These files all have the same format. For example, `mash3_datasets.csv` could look like:
```
target_sample,dataset
sampleA,sampleA
sampleA,sampleB
sampleA,sampleC
sampleB,sampleB
sampleB,sampleD
sampleB,sampleE
```
In this case, when running pipelines with `mash3`, MAG-E will supply only the samples A,B,C when binning the target sample A with the various specified binners.

- `refiners`:
`refiners` has a more complex structure. A binning refiner integrates binning outputs (each coming from potentially different pipleines) into a final output. `refiners` takes a list, each element of which has a `[REFN, REFN_OPT, [PIPELINES]]` structure, where `REFN` is the refiner name, `REFN_OPT` are the command line refiner options, and `PIPELINES` is a list of pipelines that the refiner will integrate over. For example, one may wish to evaluate DAS Tool integrated over 
CONCOCT in `mash20` mode, when run on MEGAHIT assemblies with `--no-bubble` and `--min-count 5` parameters, and METABAT2 in `single` mode, when run on default metaSPAdes assemblies. This is encoded in the element `['DAS_Tool', 'default', [['MEGAHIT', '--no-bubble --min-count 5', 'CONCOCT', 'default', 'mash20'], ['metaSPAdes', 'default', 'METABAT2', 'default', 'single']]]`. We may also wish to evaluate DAS Tool with, say, CONCOCT and METABAT2 both run in `single` mode, on default MEGAHIT. Together the complete refiners list would look like `[['DAS_Tool', 'default', [['MEGAHIT', '--no-bubble --min-count 5', 'CONCOCT', 'default', 'mash20'], ['metaSPAdes', 'default', 'METABAT2', 'default', 'single']]], ['DAS_Tool', 'default', [['MEGAHIT', 'default', 'CONCOCT', 'default', 'single'], ['MEGAHIT', 'default', 'METABAT2','default', 'single']]]]`. A bit cumbersome..., but, hopefully straightforward. 

- `qctools`:
`qctools` is a list of quality control tools that MAG-E will run on each bin. `qctools` takes a list.
For example, to run CheckM2 and GUNC, one passes `['CheckM2', 'GUNC']` to qctools.

### Verifying the project configuration.

After populating the configuration running:
```
python -m maggie verify-project
```
Checks that each file required by MAG-E is available in the config, sets up core directories,
```
TODO: LS THE CORE DIRECTORIES
```
and makes sure the assembler, binner, and refiner names specified in the profile match those
supported by the MAG-E API. This check is only performed if (`use_api` is `yes`), as
MAG-E supports evaluation of bins without using the API to allow for flexibilty in pipeline evaluation.
See "Constructing bins for evaluation" for considering whether or not to use API to construct bins.

Note how we no longer need to specify the project name or directory? MAG-E caches projects and re-loads the project config to streamline the command line. See "Project caching" below for details. 
## Required format for genome fasta
Discuss here how they need to be in a format like MGYG, so a constant genome name, an underscore, and then contig name MGYG000000512_1. Fileames should be genometoken.fasta.gz
## Simulating datasets

To simulate datasets, we firstwe make a MAG-E database of the input genomes. This involves clustering the genomes at the strain level, and constructing a Sylph sketch over them. See `--help` for options.
```
python -m maggie make-maggie-db
```

With the database constructed, we build mirror specifications of each sample, which 
lists a set of genomes from the database that closely match the sample in ANI and abundance. 
Specifications provide the input genomes for simulation and the ground truth for MAG-E evaluations.
```
python -m maggie make-mirrors
```

With the specification constructed, we can then simulate metagenomes.
```
python -m maggie simulate-mgx
```
## Constructing bins for evaluation

Now we've simulated datasets, we can construct bins using the pipelines specified in `config.yaml`. 
To do this, we first need to construct a manifest of tasks. 
```
python -m maggie construct-tasks
```

This outputs `manifest.csv` in the project directory, which enumerates all MAG-generation tasks. Each task consists of assembly, followed by binning in a particular mode, followed by potential refinement, and finally quality control. The total number of tasks is the product over this space. For example, say `config.yaml` specified three samples, two assemblers both run with two different options, two binners both run with default options in two binning modes, plus one refiner with its pipeline. That's 3 (sample) x 4 (assembly tasks) x (4 bin tasks) + 3 (sample) x 1 refiner = 51 total tasks. `manifest.csv` contains a row for each task, and each row specifies the task composition and the directories where the output needs to be written for each task in order for MAG-E to run evaluations.

### The task directory structure

Over the product space there is redundancy; all 51 tasks rely on 12 assemblies (3 samples, 4 assembly tasks), and the one refiner relies upon avaiable binning outputs. To avoid redundant computation, MAG-E caches the assemblies such that all downstream tasks can use them, and the refiner will treat the pipelines it integrates over as "cached". For tasks to find the correct cached information, MAG-E relies upon a fixed directory structure, which `construct-tasks` creates. 

Each assembly task has its own directory located in `assembly_cache`, with a subdirectory for each sample. In our example:
```
assembly_cache/
                assembly_task_1/
                                sample1/
                                sample2/
                                sample3/
                assembly_task_2/
                                ...
                assembly_task_3/
                assembly_task_4/
```
Each bin task (there are 17 in our example; 4 (assembler tasks) x 4 (bin tasks) + 1 (refiner)) has a directory in `bintask_dir` with the 
same structure:
```
bintask_dir/
            bin_task_1/
                       sample_1/
...
            bin_task_17
```

In order for MAG-E to evaluate pipelines, the output of assembly and bin tasks must be placed in the correct directories. 

### Constructing assemblies and bins

This is where the MAG-E API comes in. By default, `use_api` is `yes` (in `config.yaml`), in which case MAG-E will use the assemblers, binners, and refiners that are internally supported in its API to perform the assembly and bin tasks. MAG-E will automatically write the assembly and binning outputs to the correct locations specified in `manifest.csv`. Using the API is the most "hands-off" way to construct assemblies and bins for evaluation. First we construct assemblies with:
```
python -m maggie run-assembly
```
which populates the `assembly_cache`, and then we run the binning and refiner tasks with:
```
python -m maggie run_binning
```
After binning completes, MAG-E has written the bins for each task to the `output/bins` directory. Each bin is a separate fasta file, called `bin.<n>.fasta`.
The assemblers, binners, and refiners that MAG-E supports can be found in the `./maggie/tasks` directory of the MAG-E repository. 
Advanced users can extend the API, adding new binners, assemblers, or refiners by writing their own python classes which implement the interfaces. 

If a user decides not to use the API (i.e `use_api` is `no`; which is set with `--use_api no` when running `initialize-project`) then
the user is expected to populate the assembly_cache and bintask_dir themselves. There are two main benefits of not using the API. First,
it allows you to use any whichever compute infrastructure and workflows (e.g your own highly-parallel scripts on a slurm infrastructure) make the most sense for a large-scale assembly and binning. Second, since binning algorithm development is an active research area, it allows complete flexibility in how you decide construct bins for MAG-E evaluation. You can bin any way you want, testing any new idea you wish, so long as you: 1) come up with a name for the binning approach (e.g "`MyNewBinner`" and specifiy it in the `config.yaml` as discussed above; and, 2) you put the bins for "`MyNewBinner`" in fasta format into the `output/bins` directories that the MAG-E `manifest.csv` has issued for the "`MyNewBinner`" tasks. 

The MAG-E command `get-task-dir` returns the correct directory to put the assembly and binning output of a particular MAG-construction pipeline. It's extremely useful when not using the API. It requires a request (either `assembly` or `bin`) and a sample name (see `--help` for details). For example,
```
pth=$(python -m maggie get-task-dir assembly ERR1136644 --assembler MEGAHIT --binner CONCOCT --binning-mode all)
echo $pth
/insomnia001/depts/pmg/users/ic2465/launch/myproject/assembly_cache/assembly_task_1/ERR1136644
```
Requests the assembly task directory for the MEGAHIT assembly run with default options (`--assembler-option` not given) on sample `ERR1136644`. We returned this to the bash variable `pth` and then `echo`ed it to standard out. When requesting `assembly` the additional options are irrelevant.
The exact same command, but requesting `bin` returns the bin task directory for the MEGAHIT assembly run with default options on sample `ERR1136644` followed by CONCOCT run default in binning mode `all`:
```
pth=$(python -m maggie get-task-dir bin ERR1136644 --assembler MEGAHIT --binner CONCOCT --binning-mode all)
echo $pth
/insomnia001/depts/pmg/users/ic2465/launch/myproject/bintask_dir/bin_task_4/ERR1136644
```
This works for refiners also, with comma separating the fields within pipelines, and semi-colon separating pipelines: 
```
pth=$(python -m maggie get-task-dir bin ERR1136644 --refiner DAS_Tool --pipelines MEGAHIT,default,CONCOCT,default,all:MEGAHIT,default,METABAT2,default,single)
echo $pth
/insomnia001/depts/pmg/users/ic2465/launch/myproject/bintask_dir/bin_task_9/ERR1136644
```
`get-task-dir` makes it easy to populate the correct assembly and bin task directories when using your own scalable workflows for assembly and binning. 

### Project caching

Notice how there was no need to tell MAG-E where the project `myproject` is when it was building the project
and validating the config? This is because MAG-E holds a project cache located in `./MAG-E/.project_cache.csv`. This project cache stores each MAG-E project, and keeps track of which one is currently being worked on. We can retrieve the current project with
```
python -m maggie get-current
myproject       /Users/izaak/Library/CloudStorage/Dropbox/Documents/Projects/testingmaggie/myproject
```
which prints the current project name and its path. The caches contents looks like:
```
cat .project_cache.tsv 
name    directory       is_current
myproject       /Users/izaak/Library/CloudStorage/Dropbox/Documents/Projects/testingmaggie/myproject    True
```

Say we add another project and ask the current:
```
python -m maggie initialize-project /Users/izaak/Library/CloudStorage/Dropbox/Documents/Projects/testingmaggie anotherproject
python -m maggie get-current
anotherproject  /Users/izaak/Library/CloudStorage/Dropbox/Documents/Projects/testingmaggie/anotherproject
```
We can see the the current project has switched to `anotherproject`. Both projects are listed in the cache,
```
cat .project_cache.tsv
name    directory       is_current
myproject       /Users/izaak/Library/CloudStorage/Dropbox/Documents/Projects/testingmaggie/myproject    False
anotherproject  /Users/izaak/Library/CloudStorage/Dropbox/Documents/Projects/testingmaggie/anotherproject       True
```
with is_current set to true for `anotherproject`. 

To reset myproject to be current: 
```
python -m maggie set-current myproject
python -m maggie get-current
myproject       /Users/izaak/Library/CloudStorage/Dropbox/Documents/Projects/testingmaggie/myproject
```

This caching allows MAG-E to switch between multiple projects whilst simplifying the command-line. 

# Tutorial

First we populate the `config.yaml`.