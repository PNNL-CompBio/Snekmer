# Snekmer

Pipeline to apply encoded Kmer analysis to protein sequences

Model mode:

* **Input:** fasta protein sequences in known families
* **Output:** models of the families
  * _Evaluation output:_ assessment of model performance

Cluster mode:

* **Input:** fasta protein sequences
* **Output:** clusters of similar proteins
  * _Evaluation output:_ assessment of how well the clusters of similar proteins represent functions

## Installation

I recommend using Anaconda to create a virtual environment. Anaconda handles dependencies and versioning, which simplifies the process of installation.


### Procedure

Create a conda environment called `kmers`:

```bash
conda create -n kmers -c conda-forge -c bioconda biopython matplotlib numpy pandas seaborn snakemake scikit-learn
```

Activate the environment:

```bash
conda activate kmers
```

Install the `snekmer` package (note: git clone step is optional if you
 already have the repo cloned locally):

```bash
# clone repository if you haven't already
git clone https://github.com/PNNL-Compbio/Snekmer/Snekmer.git

# install from cloned repository
cd Snekmer
pip install .
```

The package should now be ready to use!

## Command-Line Interface

To run `snekmer`, create a `config.yaml` file containing desired
  parameters. (A template is provided at `snekmer/config.yaml`.)
  Note that the `config.yaml` file should be stored in the same
  directory as input files. Snekmer assumes that input files are
  stored in the `input` directory, and automatically creates an
  `output` directory to save all output files. Snekmer also assumes
  background files, if any, are stored in `input/background/`.
  An example of the assumed directory structure is shown below:

```
.
├── config.yaml
├── input
│   ├── background
│   │   ├── X.fasta
│   │   ├── Y.fasta
│   │   └── etc.
│   ├── A.fasta
│   ├── B.fasta
│   └── etc.
├── output
│   ├── ...
│   └── ...
```

<!-- In particular, be sure to set `output: save_dir` to the desired output file
 directory, and make sure that `input: fasta_dir` is pointing toward the
 directory containing .fasta input files. -->

### Modes

Snekmer has two operation modes: `model` (supervised modeling) and `cluster`
  (unsupervised clustering). Users may choose either mode to best suit their
  use case.

The mode must be specified in the command line, e.g.:

```bash
snekmer model [--options]
```
or

```bash
snekmer cluster [--options]
```

Once the `config.yaml` file has been updated, I recommend a dry run:

```bash
snekmer [mode] --dryrun
```

(For instance, in supervised mode, run `snekmer model --dryrun`.)

The output of the dry run shows you the files that will be created by the
 pipeline. If no files are generated, double-check   that your directory
 structure matches the format specified above.

When you are ready to process your files, run:

```bash
snekmer [mode]
```

### Partial Workflow

To execute only a part of the workflow, the parameter `--until` can be invoked.
For instance, to execute the workflow only through the kmer vector generation
step, run:

```bash
snekmer [mode] --until vectorize
```

### Extra Notes

The `snekmer` CLI is ready-to-use in the above format, but if you run
  `snekmer --help`, you'll notice many extra parameters.
  Ignore these for now; these are a WIP still!
