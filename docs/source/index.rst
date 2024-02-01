.. Snekmer documentation master file, created by
   sphinx-quickstart on Fri Jun 10 14:20:41 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Snekmer: Reduced K-mer Encoding for Protein Sequences
=====================================================

Snekmer is a software package designed to reduce the representation of protein sequences
by combining amino acid reduction (AAR) with the kmer approach. Based on the AAR-kmer representations,
Snekmer subsequently (1) clusters sequences using various unsupervised clustering algorithms,
(2) generates supervised machine learning models, or (3) searches sequences against pre-trained models
to determine probabilistic annotations. A lightweight non-machine-learning approach is also available, this 
method (1) records kmer counts from training sequences and (2) applies cosine distance calculations to new
sequences to predict the nearest annotation and generate a confidence score.  

.. image:: ../../resources/snekmer_workflow.svg
        :align: center
        :width: 700
        :alt: Snekmer workflow overview

There are 6 operation modes for Snekmer: ``cluster``, ``model``, ``search``, ``motif``, ``learn``, and ``apply``.

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

**Motif mode:** The user supplies files containing sequences in an appropriate format (e.g. FASTA). Snekmer applies the relevant workflow steps and outputs a table (.csv) for each family, which shows the SVM weight and associated p-value for each kmer.


**Learn mode:** The user supplies files containing sequences in an appropriate format (e.g. FASTA) as well as an annotation file. Snekmer generates a kmer counts matrix with the summed kmer distribution of each annotation recognized from the sequence ID. Snekmer then performs a self-evaluation to assess confidence levels. There are two outputs, a counts matrix, and a global confidence distribution. 


**Apply mode:** The user supplies files containing sequences in an appropriate format (e.g. FASTA) and
and the outputs received from Learn. Snekmer uses cosine distance to predict the annotation of each sequence from the kmer counts matrix. 
The output is a table for each file containing sequence annotation predictions with confidence levels.



.. toctree::
   :caption: Getting Started
   :maxdepth: 1
   :hidden:
   
   getting_started/install
   getting_started/config
   getting_started/cli
   getting_started/usage

.. toctree::
   :caption: Tutorial
   :maxdepth: 1
   :hidden:

   tutorial/index
   tutorial/snekmer_demo

.. toctree::
   :caption: Background
   :maxdepth: 1
   :hidden:

   background/overview
   background/ml
   background/params

.. toctree::
   :caption: Troubleshooting
   :maxdepth: 1
   :hidden:

   troubleshooting/common
