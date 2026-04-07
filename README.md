# MAG-E -- Metagenome Assembled Genome Evaluator

Welcome to MAG-E, a framework for evaluating metagenome assembled genome (MAG) construction pipelines at scale with ground truth!

This readme explains of how interact with MAG-E command-line interface, in order to:
- Construct realistic simulations of the dataset.
- Evaluate MAG-generation pipelines on the dataset.
- Evaluate MAG-generation pipelines with non-default parameters (e.g non-default assembly parameters) to explore how the parameter space impacts MAG construction performance.
- Evaluate MAG-generation pipelines on particular classes of contigs to explore how well particular genomic elements are recovered.

Interacting with MAG-E on the command-line is done by running MAG-E as a python module, like so:
```
python -m maggie --help
```
This will print the master help for MAG-E. Throughout this readme, note that each of the commands may have options which you can read about by calling `--help` on the command. 

## Installation

## Starting a new project 
MAG-E operates on "projects" which are directories with a defined structure and content. A project houses all the files required by MAG-E to simulate reads from a dataset and perform an evaluation of MAG pipelines (herein, "pipelines") on that dataset. Generally, you want a project for each dataset you wish to evaluate.

### Initializing a project
To start a new project run:
```
python -m maggie initialize-project  /absolute/path/to/launch/ myproject
```

Which constructs a new project directory, `myproject`, at `/absolute/path/to/launch/` and inside the project writes a configuration yaml file: 
```
ls ../launch/myproject 
config.yaml
```

### Configuring MAG-E - required files

The configuration file `config.yaml` is central to MAG-E's organization of the entire project. The required input files and pipelines MAG-E will evaluate are both specified through the config. After initialization, `config.yaml` looks like: 

```
project_base: /absolute/path/to/launch/myproject
project_name: myproject 
assemblers: [[ASM, OPT], [ASM, OPT]]
binners: [[BIN, OPT], [BIN, OPT]]
binning_modes: [MODE, MODE]
cluster_assignments: 'NULL'
genomes_dir: 'NULL'
prefix1: 'NULL'
qctools: ['NULL', 'NULL']
read_counts: 'NULL'
refiners: [[REFN, REFN_OPT, [[ASM, ASMOPT, BIN, BINOPT, MODE], [ASM, ASMOPT, BIN, BINOPT, MODE]]]]
samples_dir: 'NULL'
use_api: 'yes'
```

Three fields are already populated: the project base, project name, and the `use_api` flag, which set to `yes`. This flag specifies how the bins from each pipeline will be provided to MAG-E in order for it to run evaluations (See "Constructing bins for evaluation").

We must now populate the `config.yaml` with the following required inputs:
- `genomes_dir`: Absolute path to directory containing genomes that will form the MAG-E database. Each genome must be in a separate gzipped fasta file (`<genome>.fasta.gz`), where `<genome>` is the genome name. In the genome files, contig headers must have the format `><genome>_<n>`, where `n` is a unique number (e.g `>genomeA_1`, `>genomeA_2`, ...).

- `samples_dir`: Absolute path to directory containing all gzipped fastq files (prefix `.fastq.gz`). Paired end data is expected. 

- `prefix1`: Prefix for the first read file in the pair, typically either `_R1` or `_1`.

- `read_counts`: Name of a csv file specifying the number of read pairs for each sample. This must be moved inside the project directory (i.e directly under `project_base`). It should have the format:
```
sample,count
SampleA,10000
...
```

- `cluster_assignments`: The name a csv, specifying the species-level cluster assignments of each genome in `genomes_dir`. Like `read_counts` it must also be moved into the project directory. There are five required fields: `genome`, `SpeciesRepr`, `GenomeType` `N50`, and `Length`. `genome` is the name of a genome, whose file must exist at `genomes_dir/<genome>.fasta.gz`. `SpeciesRepr` lists the species-level cluster identifier for each genome. The `SpeciesRepr` values must themselves be genomes present in the `genome` columns, as they also define which genome in the species is the cluster representative. `N50` is the N50 of each genome fasta, and `Length` is the base pair length of the genome. Finally, `GenomeType` specifies whether the genome is a metagenome-assembled genome or sequenced as an isolate. It can take the values `MAG` or `Isolate`. MAG-E will only evaluate performance against `Isolate` genomes. An example of the structure for the `cluster_assignments` csv is as follows:
```
genome,SpeciesRepr,GenomeType,N50,Length
genomeA,genomeA,Isolate,25000,4000000
genomeB,genomeA,MAG,5000,1000000
genomeC,genomeA,Isolate,100000,2000000
genomeD,genomeD,Isolate,80000,2000000
```
To clarify the above descriptions, genomes A-C are in the same species cluster, with A being the representative.

### Configuring MAG-E - specifying pipelines

We now add pipelines to the configuration for MAG-E to evaluate. A pipeline consists of assembly, followed by binning in a particular mode, followed by potential refinement, and finally quality control. For each stage different algorithms that can be used, and different combinations of algorithms with their associated parameters give different pipelines. We need to specify which algorithms we wish to use for each stage. 
- `assemblers`: A list of lists. Each element of the outer-list is a 2-tuple, consisting of an assembler name followed by a string of command-line arguments used to run it. `default` is used to run the assembler with default options. For example, say we want to evaluate pipelines with MEGAHIT in default, MEGAHIT with no bubble merging but --min-count of 5, and metaSPAdes in default. We would populate `assemblers` with `[[MEGAHIT, default], [MEGAHIT, '--no-bubble --min-count 5'], [metaSPAdes, default]]`.

- `binners`: A list of lists. Follows the same format as `assemblers`. So, `[[METABAT2,default], [CONCOCT,default]]` would evaluate METABAT2 and CONCOCT MAG-pipelines in default. 

- `binning_mode`: A list of string tokens. When binning a sample (the "target sample") binning mode determines which other samples will be used. Binning mode is very general, allowing any set of user-defined samples to bin the target. Say we want to try the binning modes `single`, `all`, `mash20`, and `mash3`, where `single` uses only the target, `all` uses all samples in the dataset, `mash20` uses the 20 closest samples to the target defined by mash distance, and `mash3`, which uses the closest 3. We would first populate binning mode with `[single, all, mash20, mash']`. Then, for each token `<token>` in the list, we must place a file `<token>_datasets.csv` in the project directory. In this case the files `single_datasets.csv`, `all_datasets.csv`, `mash20_datasets.csv`, `mash3_datasets.csv`. These files must then list the samples that will be used to bin each target sample. For example, `mash3_datasets.csv` could look like:
```
target_sample,dataset
sampleA,sampleA
sampleA,sampleB
sampleA,sampleC
sampleB,sampleB
sampleB,sampleD
sampleB,sampleE
```
`target_sample` and `dataset` are the required headers, where each row specifies a target sample, and one sample that will be used to help bin the target when run in this mode. In the example, when running pipelines with mode `mash3`, MAG-E will supply only the samples A,B,C to binners when binning the target sample A.

- `refiners`: A list of lists. `refiners` has a more complex structure. A binning refiner integrates the binning outputs from different pipelines into a final output. Each element of the list passed to `refiners` has the form `[REFN, REFN_OPT, [PIPELINES]]`, where `REFN` is the refiner name, `REFN_OPT` is a string of command-line options for the refiner, and `PIPELINES` is the list of pipelines that the refiner will integrate over. For example, one may wish to run DAS Tool with default settings integrated over CONCOCT in `mash20` mode, when run on MEGAHIT assemblies with `--no-bubble` and `--min-count 5` parameters, and METABAT2 in `single` mode, when run the same assembly. This example would be encoded encoded as `[DAS_Tool, default, [[MEGAHIT, '--no-bubble --min-count 5', CONCOCT, default, mash20], [MEGAHIT, default, METABAT2, default, single]]]`. A bit cumbersome..., but, hopefully straightforward. Note that, currently, we require each pipeline to run on the same assembly, which many refiners expect. 
- `qctools`: `qctools` A list of token. This specifies the quality control tools MAG-E will run on each bin. For example, to run CheckM2 and GUNC, one passes `['CheckM2', 'GUNC']`. 

### An example configuration.
Here is an example of a complete configuration file. 
```
project_base: /insomnia001/depts/pmg/users/ic2465/launch/myproject
project_name: myproject
assemblers: [[MEGAHIT, default], [metaSPAdes, default]]
binners: [[METABAT2, default], [CONCOCT, default]]
binning_modes: [single, all]
cluster_assignments: cluster_assignments.csv
genomes_dir: /insomnia001/depts/pmg/users/ic2465/launch/myproject/genomes
prefix1: R1
qctools: [GUNC, CheckM2]
read_counts: read_counts.csv
refiners: [[DAS_Tool, default, [[MEGAHIT, default, CONCOCT, default, all], [MEGAHIT, default, METABAT2, default, single]]]]
samples_dir: /insomnia001/depts/pmg/users/ic2465/launch/myproject/small_samples
use_api: 'yes'
```

In this configuration, for all samples in `small_samples` we will run MEGAHIT and metaSPAdes in default, and bin them with METABAT2 and CONCOCT in default. MAG-E takes the product over the space of samples, binners, assemblers, plus refiners. In total there are 2 (assemblers) x 2 (binners) x 2 (binning_modes) + 1 (refiners). So nine pipelines that will be evaluated. If we add for example, `(MEGAHIT`, `--no-bubble --min-count 3)` we would now have 3 (assembler) x 2 (binners) x 2 (binning_modes) + 1 (refiners) = 13 pipelines evaluated. MAG-E therefore allows us explore how the parameter space of algorithms impact binning performance by treating the same algorithm run with different options as distinct. Pipelines with the same algorithm but different options get issues a distinct prefix (e.g, `MEGAHIT(1)`) such that they can be distinguished. MAG-E will run each pipeline over each sample. So, if we have 3 samples in the `samples_dir` we'd run a total of 27 MAG-generation tasks. 

### Verifying the project configuration.

We verify the configuration is valid by running 
```
python -m maggie verify-project
```
This checks that each file required by MAG-E is present at the correct location, sets up core directories, and makes sure the tokens used for assemblers, binners, and refiners match those supported by the MAG-E API. This token check is only performed if `use_api` is `yes`, which assumes the algorithms MAG-E currently supports in its API will be used to construct bins. See "Constructing bins for evaluation" for considering whether to use the API or not.

We no longer need to specify the project name or directory. This is because projects are cached by MAG-E in the `.project_cache.tsv` file in the MAG-E repo. MAG-E can keep track of multiple projects at once using this cache. At any time, exactly one project the cache is set to the "current project" and MAG-E loads this project's configuration each time the command-line is used to streamline commands. See "Project caching" below for details. 

## Simulating datasets

Now we have a verified configuration, we can begin simulating datasets. Running 
```
python -m maggie make-maggie-db
```
creates a MAG-E database of the genomes in `genomes_dir`. Database construction involves clustering the genomes at the strain level, and constructing Sylph sketches over the genomes. 

With the database constructed, we build mirror specifications of each sample. Each specification of a sample lists a set of genomes from the database that closely match it in ANI and abundance. Specifications define the set of genomes and their abundance that will be used to simulate samples, and provide ground truth for MAG-E evaluations.
```
python -m maggie make-mirrors
```

With the specification constructed, we can then simulate metagenomes. 
```
python -m maggie simulate-mgx
```
For each of these commands, check their `--help` to see possible options, for example, increasing the thread count. 

## Constructing bins for evaluation

Now we've simulated datasets we can start running pipelines and constructing MAGs. There are potentially many pipelines, and many MAG-generation tasks that MAG-E needs to keep track of in order to successfully perform evaluations. To keep track of all this information, MAG-E constructs a manifest, each row of which describes a particular MAG-generation task: which assembler, binner, refiner, and binning mode is running on which sample, and in which directories of the project the associated assemblies, binning outputs, quality control output, and MAG-E raw metrics can be found. We construct this manifest with

```
python -m maggie construct-tasks
```

which outputs `manifest.csv` in the project directory.

### The MAG-generation task directory structure
Running `construct-tasks` also constructs the directories where MAG-E will put the various outputs of the MAG-generation tasks. These directories have a defined structure, which MAG-E relies upon in order to correctly locate files. To avoid redundant computation, a MAG-generation task is decomposed into two stages: the assembly task and the binning task portions. This separation allows MAG-generation tasks that share same assembler (and assembler option) to use the same assembly files, rather than recomputing them every time. Each assembly task has its own directory located in `assembly_cache`, with a subdirectory for each sample. For the example `config.yaml` above where we run MEGAHIT and metaSPAdes in default, assuming three samples called `sampleA`, `sampleB` and `sampleC`, `assembly_cache` contains:
```
assembly_cache/
    assembly_task_1/
                sampleA/
                sampleB/
                sampleC/
    assembly_task_2/
                sampleA/
                sampleB/
                sampleC/
```
Each of the leaf directories e.g `assembly_cache/assembly_task_1/sampleA`, should contain the assembly.

The bin tasks rely upon the assemblies cached in `assembly_cache`. The number of bin tasks equals the number of MAG-generation tasks. In our example, there are nine total MAG-generation tasks, and so, there are nine bin task directories in `bintask_dir`:
```
bintask_dir/
    bin_task_1/
        sample_A/
        sample_B/
        sample_C/
    bin_task_2/
...
```

### Constructing assemblies and bins

In order for MAG-E to perform evaluations of the various pipelines, we now need to output the assemblies and bins from each pipeline into the correct directories such that MAG-E can keep track. The `manifest.csv` already tells us which `assembly_cache` and `bintask_dir` subdirectory to place each assembly and binning output from every MAG-generation task for MAG-E to work successfully. There are two ways to use it to construct assemblies and bins. 

The first is to use the MAG-E API, which is used when `use_api` is set to `yes` in the `config.yaml`. By default the API is used. In this case, MAG-E will use the internally supported assemblers, binners, and refiners in order to run the MAG-generation tasks and automatically populate the `assembly_cache` and `bintask_dir` according to the manifest. The major benefit is that this is the most "hands-off" way to construct assemblies, bins, and correctly populate the directories. The downside is that, currently, each MAG-generation task is run sequentially, and so depending on your compute infrastructure this can be slow if you have many tasks. The supported algorithms can be found in the repo subdirectory `./maggie/tasks`, and are simple python class interfaces that wrap each algorithm, run it on the command-line, and format the outputs to work with the downstream part of MAG-E. These were built with the idea that it would be straightforward expand the automated part of assembly and bin construction with ease (see Developer Section to understand the interface and how to extend). 

Using the API is straightforward. To run the assembly tasks:
```
python -m maggie run-assembly
```
which populates the `assembly_cache` according to the manifest. Many binners require read-to-contig mapping (or count) information as input to their algorithms. We compute the `bam` files with 
```
python -m maggie run-mapping
```
which maps every `dataset` sample defined in a binning mode csv against every `target_sample` and places these files in the `assembly_cache`. 

Then, to run the binning and refiner tasks:
```
python -m maggie run_binning
```
which populates the `bintask_dir` according to the manifest. 

The second way to use the manifest to construct assemblies and bins is to query it to tell you the correctly location to put assemblies and bins, but put them there yourself, without using the API. 
There a two main benefits for not using the API (`use_api` is `no`). First, it provides the most flexibility in how you construct the assemblies and bins. Many users have powerful compute infrastructures available to them, but with different architectures, and accepted workflows. Not using the API allows you to make use of whichever workflows and compute infrastructure makes the most sense for you to efficiently run lots of assemblies and binning tasks (e.g your own in-house slurm or sun grid workflows, aws batch jobs, or nextflow or snakemake scripts). The second major and very powerful advantage is that, you can bin in any way you want, testing any new idea you wish, and evaluate it with MAG-E. Since binning research is an active field, the API interface we offer may not fit a new approach. So long as you come up with a new `binner` token for the approach, e.g `MyCrazyBinner`, and add it to the config as discussed above, MAG-E can keep track of it, make `assembly_cache` and `bintask_dir` directories for it, and ultimately evaluate it.

The main downside is that you have to put the assembly and binning outputs in the correct locations, in the correct format. To help with this MAG-E implements the command `get-task-dir` which queries the the manifest and returns the correct directory to write an assembly or binning output for a particular MAG-generation task. It's extremely useful when not using the API for putting stuff in the correct place. It requires a request (either `assembly` or `bin`) and a sample name (see `--help` for details). Running with the example `config.yaml` from above. Suppose we want the correct assembly task directory for the MAG-generation task which runs `MEGAHIT` in default, `CONCOCT` in default, and binning mode `all` on sample `ERR113664`:
```
pth=$(python -m maggie get-task-dir assembly ERR1136644 --assembler MEGAHIT --binner CONCOCT --binning-mode all)
echo $pth
/insomnia001/depts/pmg/users/ic2465/launch/myproject/assembly_cache/assembly_task_1/ERR1136644
```
Note that, since this is requesting the `assembly` task, one could drop the binner, and binning mode and get the same return value. Also, since `MEGAHIT` is run in default, the `--assembler-option` part of the query is not needed. The exact same query, but instead requesting `bin`, returns the bin task directory for MAG-generation task:
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

For each assembly task, the contigs need to be added to the correct `assembly_cache` subdirectory as a fasta file with filename format, `<sample>.fasta`, e.g `ERR1136644.fasta` above, and each contig should be called `>NODE_<n>` where `<n>` is a unique numerical value. For each bin task, each bin should be written to the correct `bintask_dir` subdirectory as a separate fasta file with filename format `bin.<n>.fasta` where `<n>` is a unique numerical value (they don't have to be consecutive). Bam files needed by binning algorithms are put into the assembly cache. For each unique pair of (target_sample, dataset sample) across all the binning mode csvs, a coordinate-sorted bam file `<target>_<dataset_sample>.bam` needs to be added to the correct `assembly_cache` subdirectory. For example `ERR1136644_ERR1136644.bam`, which maps the reads of sample `ERR1136644` to its assembly would be needed for single-sample binning of the sample. 

## Running quality control
We can now run quality control algorithms on bins to predict their quality with:
```
python -m maggie run-quality-control
``` 

Currently MAG-E supports CheckM2 and GUNC. This follows a similar pattern to the assembly, binner, and refiner API (see Developer Section). Currently, MAG-E does not support non-API quality control.

## MAG-pipeline evaluation

We now have nearly everything we need to evaluate pipelines with MAG-E. One final file set we need to compute is the ground truth assignments of every contigs to every genome. We do this with:
```
python -m maggie construct-ground-truth
```
This writes `<sample>_gt_table.csv` files to the correct locations in `assembly_cache`. Each file maps the contigs to their genome(s) of origin, providing ground truth for MAG-E evaluation. 

We can now calculate for each MAG-generation task, the ground truth metrics (precision, and recall, and F-score) for each genome. We do this with:
```
python -m maggie calc-per-genome-metrics
``` 

This will generate a number of files, the most important of which are prefixed with `.genome_metrics.parquet` and will be written to the `evaluation` directory. This binary files report the ground truth metrics of each genome from each MAG-generation task. There's one file per (assembler,binner) combination, to avoid the files becoming overly large. 

To get the final MAG-E evaluation, we call:
```
python -m maggie evaluate-pipelines
```
which first constructs the recoverable genome set (recall >= 0.7, precision >= 0.9 by default), and then builds a linear-mixed model of MAG-pipeline performance using the metrics from the recoverable set. MAG-E reports the estimated marginal means for each pipeline, and various performance comparisons between assemblers, binners and refiners, and binning modes to a set of files in `evaluation/LMM`:
- `all_pipelines_fscore.csv`: Reports the estimated marginal mean F-score for each MAG-pipeline.
- `all_pipelines_precision.csv`: Reports the estimated marginal mean (1-precision) for each MAG-pipeline.
- `all_pipelines_recall.csv`: Reports the estimated marginal mean recall for each MAG-pipeline.
- `pipeline_categories*.csv`: Report the estimated marginal mean of particular pipeline categories, e.g the mean of all pipelines running `MEGAHIT`. 
- `all_pairs*.csv`: Statistical difference tests between all pairs of pipelines for each metric. 

Note that tests over categories are currently printed to stdout. 
With these tables and printed outputs, one can plot and compare the different MAG-pipelines. 

## Contig level evaluations with MAG-E

With MAG-E we can evaluate MAG-pipeline performance on particular classes of contigs. Currently MAG-E only supports this for pipelines run with default assembly and binning command-line options. To evaluate contig classes for each assembly, you must construct a "contig properties" file. These files have the following format:
```
sample,contig,prop_d_<x1>,prop_d_<x2>,...
```
There must be one row per contig. The two required columns `sample`and `contig` must specify the name of the sample and contig that each row represents. Contig properties must have the name format `prop_d_<xn>` where `<xn>` is the name of the property. For example, it could be `prop_d_prophage` which has two levels, `prophage` and `not_prophage`, which is used to classify each contig as containing a prophage annotation or not. This will allow MAG-E to evaluate the pipeline performance on contigs with a prophage annotation. For example, say `assembly_cache/assembly_task_1/ERR1136644` is the correct directory for running metaSPAdes in default on sample `ERR1136644`. In this directory, we need to add a file with filename format `<sample>_contig_properties.csv`, i.e `ERR1136644_contig_properties.fasta`, which may look like
```
sample,contig,prop_d_prophage,prop_d_containsOri,prop_d_NumberOfPolyARepeats
ERR1136644,NODE_1,no_prophage,no_ori,0
ERR1136644,NODE_2,no_prophage,no_ori,1
ERR1136644,NODE_3,prophage,no_ori,2
ERR1136644,NODE_4,no_prophage,ori,0
ERR1136644,NODE_5,no_prophage,ori,1
```
Which could, for example, check the performance of contigs with prophage annotations, 
origin of replication annotations, or with different numbers of A-repeats of length > 10, 
however the user should define it.

## Project caching
MAG-E caches projects, allowing the user to work on multiple projects with a streamlined command line, and switch between them. MAG-E does this by caching project names and filesystem locations in `./MAG-E/.project_cache.csv`. One project has a "current" status, which MAG-E will automatically load the configuration of when running commands. We can retrieve the current project with
```
python -m maggie get-current
myproject       /absolute/path/to/launch/myproject
```
which prints the current project name and its path. The cache contents looks like:
```
cat .project_cache.tsv 
name    directory       is_current
myproject       /absolute/path/to/launch/myproject    True
```

Say we add another project and ask the current:
```
python -m maggie initialize-project /absolute/path/to/launch/ anotherproject
python -m maggie get-current
anotherproject  /absolute/path/to/launch/anotherproject
```
We see that the current project has switched to `anotherproject`. Both projects are listed in the cache,
```
cat .project_cache.tsv
name    directory       is_current
myproject       /absolute/path/to/launch/myproject    False
anotherproject  /absolute/path/to/launch/anotherproject       True
```
with is_current set to true for `anotherproject`. 

To reset `myproject` to be current: 
```
python -m maggie set-current myproject
python -m maggie get-current
myproject       /absolute/path/to/launch/myproject
```
## Developer section
To be added. 
