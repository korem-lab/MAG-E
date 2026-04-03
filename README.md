# MAG-E -- Metagenome Assembled Genome Evaluator

Welcome to MAG-E, a framework for evaluating metagenome assembled genome (MAG) construction pipelines at scale with ground truth!

This readme provides an explaination of how to run MAG-E on real metagenomic datasets to: 
- Construct realistic simulations of the dataset
- Evaluate MAG-generation pipelines on the simulated dataset dataset
## Installation

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
It has four required fields: `genome`, `SpeciesRepr`, and `GenomeType` and `N50`. `genome` is the name 
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
`refiners` has a more complex structure. A binning refiner integrates binning outputs (each coming from potentially different pipleines) into a final output. `refiners` takes a list, each element of which has a `[REFN, [PIPELINES]]` structure, where `REFN` is the refiner name, and `PIPELINES` is a list of pipelines that the refiner will integrate over. For example, one may wish to evaluate DAS Tool integrated over 
CONCOCT in `mash20` mode, when run on MEGAHIT assemblies with `--no-bubble` and `--min-count 5` parameters, and METABAT2 in `single` mode, when run on default metaSPAdes assemblies. This is encoded in the element `['DAS_Tool', [['MEGAHIT', '--no-bubble --min-count 5', 'CONCOCT', 'default', 'mash20'], ['metaSPAdes', 'default', 'METABAT2', 'default', 'single']]]`. We may also wish to evaluate DAS Tool with, say, CONCOCT and METABAT2 both run in `single` mode, on default MEGAHIT. Together the complete refiners list would look like `[['DAS_Tool', [['MEGAHIT', '--no-bubble --min-count 5', 'CONCOCT', 'default', 'mash20'], ['metaSPAdes', 'default', 'METABAT2', 'default', 'single']]], ['DAS_Tool',[['MEGAHIT', 'default', 'CONCOCT', 'default', 'single'], ['MEGAHIT', 'default', 'METABAT2','default', 'single']]]]`. A bit cumbersome..., but, hopefully straightforward. 

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

### Simulating datasets

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

### Constructing bins for evaluation

MAG-E can bin samples using the assemblers, binners, and refiners that are internally supported by its
API. For flexibility, it can also evaluate bins constructed by means not supported by the API. See "Constructing bins for evaluation" section for the considerations over whether to use the API. 
If the API is used (`use_api` is `yes`), the verification checks whether the assemblers, binners, and refiners specified in `config.yaml` are supported. To get the list of currently supported tools, run `python -m maggie list-supported`. If `use_api` is `no`, this check is not run. In either case, MAG-E will check whether the quality control tools specified in `config.yaml` are supported. 



# Tutorial

First we populate the `config.yaml`.