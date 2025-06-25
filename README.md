# Snekmer: A scalable pipeline for protein sequence fingerprinting using amino acid recoding (AAR)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7662597.svg)](https://doi.org/10.5281/zenodo.7662597)
[![CI](https://github.com/PNNL-CompBio/Snekmer/actions/workflows/action.yml/badge.svg)](https://github.com/PNNL-CompBio/Snekmer/actions)
[![Documentation Status](https://readthedocs.org/projects/snekmer/badge/?version=latest)](https://snekmer.readthedocs.io/en/latest/?badge=latest)
[![Snakemake](https://img.shields.io/badge/snakemake-=7.0.0-brightgreen.svg?style=flat)](https://snakemake.readthedocs.io)

Snekmer is a software package designed to reduce the representation of protein sequences
by combining amino acid reduction (AAR) with the kmer approach. Based on the AAR-kmer representations,
Snekmer subsequently (1) clusters sequences using various unsupervised clustering algorithms,
(2) generates supervised machine learning models, or (3) searches sequences against pre-trained models
to determine probabilistic annotations.

<p align="center">
  <img align="center" src="resources/snekmer_workflow.svg">
</p>

There are six operation modes for Snekmer: `cluster`, `model`, and `search`, `learn`, `apply`, and `motif`.

**Cluster mode:** The user supplies files containing sequences in an appropriate format (e.g. FASTA).
Snekmer applies the relevant workflow steps and outputs the resulting clustering results in tabular form (.CSV),
as well as the cluster object itself (.cluster). Figures are also generated (e.g. t-SNE, UMAP) to help the user
contextualize their results.

**Model mode:** The user supplies files containing sequences in an appropriate format (e.g. FASTA).
Snekmer applies the relevant workflow steps and outputs the resulting models as objects (.model). Snekmer also
displays K-fold cross validation results in the form of figures (AUC ROC and PR AUC curves) and a table (.CSV).

**Search mode:** The user supplies files containing sequences in an appropriate format (e.g. FASTA)
and the models they wish to search their sequences against. Snekmer applies the relevant workflow steps
and outputs a table for each file containing model annotation probabilities for the given sequences.

**Learn mode:** The user supplies files containing sequences in an appropriate format (e.g. FASTA)
as well as an annotation file. Snekmer generates a kmer counts matrix with the summed kmer distribution
of each annotation recognized from the sequence ID. Snekmer then performs a self-evaluation to assess
confidence levels. There are two outputs, a counts matrix, and a global confidence distribution.

**Apply mode:** The user supplies files containing sequences in an appropriate format (e.g. FASTA)
and the outputs received from Learn. Snekmer uses cosine distance to predict the annotation of each
sequence from the kmer counts matrix. The output is a table for each file containing sequence annotation
predictions with confidence levels.

**Motif mode:** The user supplies files containing sequences in an appropriate format (e.g. FASTA)
and the outputs received from Model. Snekmer performs a feature selection workflow to produce a 
list of motifs ordered by degree of conservation and a classification model using the selected features (.model).

## How to Use Snekmer

For installation instructions, documentation, and more, refer to
the [official documentation](https://snekmer.readthedocs.io).

To run the demonstration example, see 
[resources/tutorial/demo_example](https://github.com/PNNL-CompBio/Snekmer/tree/main/resources/tutorial/demo_example).


Snekmer was written and is maintained by the following PNNL development team: Christine Chang, Jeremy Jacobson, Abby Jerger, Bill Nelson, and Jason McDermott.

## Citation Guidance

1. McDermott, Jason E., Chang, Christine H., Jerger, Abby, Nelson, William B., & Jacobson, Jeremy R. (2023). Snekmer: A scalable pipeline for protein sequence fingerprinting using amino acid recoding (AAR) (v1.0.3). Zenodo. [https://doi.org/10.5281/zenodo.7662597](https://doi.org/10.5281/zenodo.7662597)
2. Christine H Chang, William C Nelson, Abby Jerger, Aaron T Wright, Robert G Egbert, Jason E McDermott, Snekmer: a scalable pipeline for protein sequence fingerprinting based on amino acid recoding, Bioinformatics Advances, Volume 3, Issue 1, 2023, vbad005, [https://doi.org/10.1093/bioadv/vbad005](https://doi.org/10.1093/bioadv/vbad005).

## License

    This material was prepared as an account of work sponsored by an agency of the United States Government. Neither the United States Government nor the United States Department of Energy, nor Battelle, nor any of their employees, nor any jurisdiction or organization that has cooperated in the development of these materials, makes any warranty, express or implied, or assumes any legal liability or responsibility for the accuracy, completeness, or usefulness or any information, apparatus, product, software, or process disclosed, or represents that its use would not infringe privately owned rights.

    Reference herein to any specific commercial product, process, or service by trade name, trademark, manufacturer, or otherwise does not necessarily constitute or imply its endorsement, recommendation, or favoring by the United States Government or any agency thereof, or Battelle Memorial Institute. The views and opinions of authors expressed herein do not necessarily state or reflect those of the United States Government or any agency thereof.

    PACIFIC NORTHWEST NATIONAL LABORATORY operated by BATTELLE for the UNITED STATES DEPARTMENT OF ENERGY under Contract DE-AC05-76RL01830
