Modes
=====

Snekmer has three operation modes: ``model`` (supervised modeling),
``cluster`` (unsupervised clustering), and ``search`` (application
of model to new sequences). We will call first two **learning modes**
due to their utility in learning relationships between protein family
input files. Users may choose a mode to best suit their use case.

The mode must be specified in the command line, e.g. to specify the
``model`` mode, the following should be called:

.. code-block:: bash

    snekmer model [--options]

In the **resources** directory, two example configuration files are included:

  - **resources/config.yaml**: Configuration file for ``snekmer model`` and ``snekmer cluster`` modes.
  - **resources/search.yaml**: Configuration file for ``snekmer search`` mode. Note that the Snekmer CLI automatically assumes that the configuration file will be named *config.yaml*, so to use the provided file, use ``snekmer search --configfile search.yaml``

.. code-block:: bash

    snekmer [mode] --dryrun

(For instance, in supervised mode, run ``snekmer model --dryrun``.)

The output of the dry run shows you the files that will be created by the
pipeline. If no files are generated, double-check   that your directory
structure matches the format specified above.

When you are ready to process your files, run:

.. code-block:: bash

    snekmer [mode]

Output
------

Each step in the Snekmer pipeline generates its own associated output files.
Both operation modes will preprocess parameters, generate labels, and
vectorize sequences based on labels. The associated output files can be
found in the respective directories.

The following output directories and files are created in both operation modes:

.. code-block:: console

    .
    ├── input/
    │   ├── A.fasta
    │   └── B.fasta
    └── output/
        ├── processed/
        │   ├── A.json             # processed parameter values for A
        │   ├── B.json             # processed parameter values for B
        │   ├── A_description.csv  # summary of sequences in A.fasta
        │   └── B_description.csv  # summary of sequences in B.fasta
        ├── labels/
        │   ├── A.txt              # kmer labels for A
        │   └── B.txt              # kmer labels for B
        ├── features/
        └── ...

Model Mode
----------

Executing ``snekmer model`` produces the following output files
and directories in addition to the files described previously.

.. code-block:: console

    .
    └── output/
        ├── ...
        ├── features/
        │   ├── A/            # kmer vectors in A kmer space
        │   │   ├── A.json.gz
        │   │   └── B.json.gz
        │   └── B/            # kmer vectors in B kmer space
        │       ├── A.json.gz
        │       └── B.json.gz
        ├── score/
        │   ├── A.pkl         # A sequences, scored
        │   ├── B.pkl         # B sequences, scored
        │   └── weights/
        │       ├── A.csv.gz  # kmer score weights in A kmer space
        │       └── B.csv.gz  # kmer score weights in B kmer space
        └── model/
            ├── A.pkl         # (A/not A) classification model
            ├── B.pkl         # (B/not B) classification model
            ├── results/      # cross-validation results table
            │   ├── A.csv
            │   └── B.csv
            └── figures/      # cross-validation results figures
                ├── A/
                └── B/

Cluster Mode
------------

Executing ``snekmer cluster`` produces the following output files
and directories in addition to the files described previously.

.. code-block:: console

    .
    └── output/
        ├── ...
        ├── features/
        │   └── full/     # kmer vectors in full kmer space for (alphabet, k)
        │       ├── A.json.gz
        │       └── B.json.gz
        └── cluster/
            ├── A.pkl     # A cluster model
            ├── B.pkl     # B cluster model
            └── figures/  # cluster figures (t-SNE)
                ├── A/
                └── B/

Search Mode
-----------

The ``snekmer search`` mode assumes that the user has pre-generated
family models using the `snekmer model` workflow, and thus operates
as an independent workflow. The location of the basis sets, scorers,
and models must be specified in the configuration file (see:
**resources/search.yaml**).

For instance, say that the above output examples have already been
produced. The user would then like to search a set of unknown
sequences against the above families.

In a separate directory, the user should place files in an input
directory with the appropriate YAML file. The assumed input file
structure is as follows:

.. code-block:: console

    .
    ├── search.yaml
    ├── input/
    │   ├── unknown_1.fasta
    │   ├── unknown_2.fasta
    │   └── etc.
    ├── output/
    │   ├── ...
    │   └── ...

The user should then modify **search.yaml** to point toward the
appropriate basis set, scorer, and model directories.

Executing ``snekmer search --configfile search.yaml`` produces the
following output files and directories in addition to the files
described previously.

.. code-block:: console

    .
    └── output/
        ├── features/
        │   ├── A/
        │   │   ├── unknown_1.json.gz
        │   │   └── unknown_2.json.gz
        │   └── B/
        │       ├── unknown_1.json.gz
        │       └── unknown_2.json.gz
        └── search/
            ├── A.csv  # A probabilities and predictions for unknown sequences
            └── B.csv  # B probabilities and predictions for unknown sequences