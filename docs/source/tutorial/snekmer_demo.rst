.. _model-cluster-search-tutorial:

Model / Cluster / Search Tutorial
==================================

This tutorial walks through three Snekmer modes — **Cluster**, **Model**, and **Search** —
using the demo data included in the repository.

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Mode
     - What it does
   * - ``cluster``
     - Unsupervised clustering of sequences by k-mer profile. Produces a cluster assignment
       table and optional figures (t-SNE, UMAP, PCA).
   * - ``model``
     - Trains one-vs-rest supervised models from annotated sequences. Produces per-family
       classification models and K-fold cross-validation figures (AUC ROC, PR AUC).
   * - ``search``
     - Scores new sequences against models from ``snekmer model``. Produces per-family
       annotation probability tables.


Demo data
---------

Input files are in ``resources/demo_sequences/model_cluster_search_inputs/``:

.. code-block:: text

   model_cluster_search_inputs/
   ├── nirS.faa        ← nitrite reductase family sequences
   ├── nxrA.faa        ← nitrite oxidoreductase family sequences
   └── TIGR03149.faa   ← TIGRFAM 03149 family sequences

All commands below assume you are running from the demo workspace directory:

.. code-block:: bash

   cd resources/model_cluster_search_demo


Configuration
-------------

A ``config.yaml`` is required in the working directory. The demo includes one
pre-configured for these three input families:

.. code-block:: bash

   cat config.yaml

Key parameters:

.. list-table::
   :header-rows: 1
   :widths: 25 20 55

   * - Parameter
     - Default
     - Description
   * - ``k``
     - ``8``
     - K-mer length
   * - ``alphabet``
     - ``2`` (solvacc)
     - Amino acid reduction alphabet (0–5 or name)
   * - ``cluster.method``
     - ``agglomerative-jaccard``
     - Clustering algorithm
   * - ``cluster.params.distance_threshold``
     - ``0.92``
     - Jaccard distance threshold for cluster splitting
   * - ``model.cv``
     - ``5``
     - K-fold cross-validation folds

See :doc:`../getting_started/config` for the full parameter reference.


Running the demo
----------------

The ``run_demo.py`` script resets the workspace, copies inputs, and runs all three modes
in order:

.. code-block:: bash

   python run_demo.py

Or run each mode individually after copying the input files and config manually:

Step 1: Copy inputs and config
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   cp ../demo_sequences/model_cluster_search_inputs/nirS.faa       input/
   cp ../demo_sequences/model_cluster_search_inputs/nxrA.faa        input/
   cp ../demo_sequences/model_cluster_search_inputs/TIGR03149.faa   input/
   cp ../config.yaml ./config.yaml

Step 2: Cluster
~~~~~~~~~~~~~~~~

.. code-block:: bash

   snekmer cluster --scheduler greedy --configfile=./config.yaml

Step 3: Model
~~~~~~~~~~~~~~

.. code-block:: bash

   snekmer model --scheduler greedy --configfile=./config.yaml

Step 4: Collect model artifacts (required before Search)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``snekmer search`` needs the per-family ``.model``, ``.kmers``, and ``.scorer``
files collected into a single directory (``output/example-model/`` by default):

.. code-block:: bash

   mkdir -p output/example-model
   cp output/model/*.model     output/example-model/
   cp output/kmerize/*.kmers   output/example-model/
   cp output/scoring/*.scorer  output/example-model/

Step 5: Search
~~~~~~~~~~~~~~~

.. code-block:: bash

   snekmer search --scheduler greedy --configfile=./config.yaml


Output structure
----------------

After running all three modes the output directory contains:

.. code-block:: text

   output/
   ├── cluster/
   │   ├── snekmer.csv               ← cluster assignment table (one row per sequence)
   │   └── figures/                  ← t-SNE / UMAP / PCA plots (if enabled)
   ├── kmerize/
   │   ├── nirS.kmers
   │   ├── nxrA.kmers
   │   └── TIGR03149.kmers
   ├── model/
   │   ├── nirS.model
   │   ├── nxrA.model
   │   └── TIGR03149.model
   ├── scoring/
   │   ├── sequences/                ← per-family sequence probability scores (.csv.gz)
   │   ├── weights/                  ← per-family k-mer weights (.csv.gz)
   │   └── *.scorer                  ← scorer objects
   ├── search/                       ← search results (one CSV per family)
   ├── Snekmer_Cluster_Report.html
   └── Snekmer_Model_Report.html


Reading cluster results
-----------------------

``output/cluster/snekmer.csv`` has one row per input sequence:

.. code-block:: python

   import pandas as pd
   df = pd.read_csv("output/cluster/snekmer.csv")
   print(df.head())

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Column
     - Description
   * - ``sequence_id``
     - Sequence identifier from the FASTA header
   * - ``filename``
     - Source FASTA file (used as family label)
   * - ``cluster``
     - Integer cluster assignment
   * - ``kmer_*``
     - K-mer feature values (one column per k-mer)


Reading model / search results
-------------------------------

For each family, ``output/scoring/sequences/<family>.csv.gz`` contains
the probability score for every sequence:

.. code-block:: python

   import pandas as pd

   family = "nirS"
   df = pd.read_csv(f"output/scoring/sequences/{family}.csv.gz")
   print(df[["sequence_id", "label", f"{family}_score"]].head())

For search results, ``output/search/<family>/<family>.csv`` contains
the annotation probability for unknown sequences:

.. code-block:: python

   family = "nirS"
   df = pd.read_csv(f"output/search/{family}/{family}.csv")
   print(df.head())

HTML summary reports for Cluster and Model are written to the output root
(``Snekmer_Cluster_Report.html``, ``Snekmer_Model_Report.html``).


Jupyter notebook
----------------

An interactive version of this tutorial (with result visualisations) is available at
``docs/source/tutorial/snekmer_demo_notebook.ipynb``.
