# KmerPipeline
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
2. Parse annotations from protein sequences into usable form (if neessary)
3. Run pipeline above with some parameters
4. Assess how sensitive the method is - that is how likely is it that a members of a cluster have the annotation that is the most prevalent for that cluster
5. Assess the specificity - that is how likely is it that an annotation maps to a single cluster

## Installation

I recommend using Anaconda to create a virtual environment. Anaconda handles dependencies and versioning, which simplifies the process of installation.


### Procedure

Create a conda environment called `kmers`:

```bash
conda create -n kmers -c conda-forge -c bioconda biopython numpy pandas snakemake scikit-learn
```

Activate the environment:

```bash
conda activate kmers
```

Install the `kmerfeatures` package (note: git clone step is optional if you
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

To run `kmerfeatures`, make sure to modify `kmerfeatures/config.yaml` to set
 the desired parameters for analysis. Then, navigate to the directory
 containing `config.yaml`. If you are using the default config file, this entails:

 ```bash
cd kmerfeatures
 ```

In particular, be sure to set `output: save_dir` to the desired output file
 directory, and make sure that `input: fasta_dir` is pointing toward the
 directory containing .fasta input files.

Once the config file has been updated, I recommend running the following:

```bash
kmerfeatures --dryrun
```
The output of the dry run shows you the files that will be created by the
 pipeline. If no files are generated, double-check `input: fasta_dir` in the
 config file, and make sure that the desired outputs are being generated.

When you are ready to process your files, run:

```bash
kmerfeatures --cores 1
```

### Extra Notes

The `kmerfeatures` CLI is ready-to-use in the above format, but if you run `kmerfeatures --help`, you'll notice many extra parameters. Ignore these for now; these are a WIP still!
