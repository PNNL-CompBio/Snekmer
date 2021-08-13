# Snekmer
Pipeline to apply encoded Kmer analysis to protein sequences

* **Input:** fasta protein sequences
* **Output:** clusters of similar proteins
  * _Evaluation output:_ assessment of how well the clusters of similar proteins represent functions

#### Steps envisioned:

1. Kmer feature generation (using KmerFeatures.py)
2. Similarity analysis - based on kmer vectors - Jaccard similarity?
3. Clustering of similarity graph to communities - MCL? Other?

#### Evaluation pipeline:
1. Start with well-annotated protein sequences
2. Parse annotations from protein sequences into usable form (if necessary)
3. Run pipeline above with some parameters
4. Assess how sensitive the method is - that is how likely is it that a members of a cluster have the annotation that is the most prevalent for that cluster
5. Assess the specificity - that is how likely is it that an annotation maps to a single cluster

## Installation

I recommend using Anaconda to create a virtual environment. Anaconda handles dependencies and versioning, which simplifies the process of installation.


### Procedure

Create a conda environment called `kmers`:

```bash
conda create -n kmers -c conda-forge -c bioconda biopython matplotlib numpy pandas snakemake scikit-learn
```

Activate the environment:

```bash
conda activate kmers
```

Install the `snekmer` package (note: git clone step is optional if you
 already have the repo cloned locally):

```bash
# clone repository if you haven't already
git clone https://github.com/biodataganache/KmerPipeline.git

# install from cloned repository
cd KmerPipeline
git checkout christine
pip install .
```

The package should now be ready to use!

## Command-Line Interface

To run `snekmer`, create a `config.yaml` file containing desired
  parameters. (A template is provided at `snekmer/config.yaml`.)
  Note that the `config.yaml` file should be stored in the same directory
  as input files. Snekmer assumes that input files are stored in the `input`
  directory, and automatically creates an `output` directory to save all output
  files. An example of the assumed directory structure is shown below:

```
.
├── config.yaml
├── input
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

Once the `config.yaml` file has been updated, I recommend a dry run:

```bash
snekmer --dryrun
```
The output of the dry run shows you the files that will be created by the
 pipeline. If no files are generated, double-check that your directory
 structure matches the format specified above.

When you are ready to process your files, run:

```bash
snekmer
```

### Partial Workflow

To execute only a part of the workflow, the parameter `--until` can be invoked.
For instance, to execute the workflow only through the kmer vector generation
step, run:

```bash
snekmer --until standardize_kmers
```

### Extra Notes

The `snekmer` CLI is ready-to-use in the above format, but if you run
  `snekmer --help`, you'll notice many extra parameters.
  Ignore these for now; these are a WIP still!
