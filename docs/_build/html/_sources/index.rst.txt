.. Snekmer documentation master file, created by
   sphinx-quickstart on Fri Jun 10 14:20:41 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Snekmer: K-mer Encoding for Protein Sequences
=============================================

Snekmer is a pipeline to apply encoded k-mer analysis to protein sequences for machine learning.

.. image:: ../../resources/snekmer_workflow.png
        :align: center
        :width: 700
        :alt: Snekmer workflow overview

There are 3 operation modes for Snekmer: ``model``, ``cluster``, and ``search``.

Model mode:
-----------

* **Input:** FASTA files containing protein sequences in known families
* **Output:** Models of the known protein families based on kmer vector analysis
   - *Evaluation output:* Assessment of model performance

Cluster mode:
-------------
* **Input:** FASTA files containing protein sequences
* **Output:** Clusters of similar proteins
   - *Evaluation output:* Assessment of how well the clusters of similar proteins represent functions

Search mode:
------------
* **Input:** FASTA files containing protein sequences; Trained model (output from snekmer model mode)
* **Output:** Predictions of family membership

.. toctree::
   :caption: Getting Started
   :maxdepth: 1
   :hidden:
   
   getting_started/install
   getting_started/cli
   getting_started/modes

.. toctree::
   :caption: Tutorial
   :maxdepth: 1
   :hidden:

   tutorial/index