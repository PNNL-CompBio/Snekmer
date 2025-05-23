Command Line Interface
======================

To run any of the five Snekmer operation modes, simply call:

.. code-block:: bash

    snekmer {mode}

Each mode has its own mode-specific options and parameters to be specified
on the command line or the ``config.yaml`` file, respectively.

For an overview of Snekmer usage, reference the help command (``snekmer --help``).

.. code-block:: console

    $ snekmer --help
    usage: snekmer [-h] [-v] {cluster,model,search,learn,apply,motif} ...

    Snekmer: A tool for kmer-based sequence analysis using amino acid reduction (AAR)

    options:
    -h, --help            show this help message and exit
    -v, --version         print version and exit

    mode:
    Snekmer mode

    {cluster,model,search,learn,apply,motif}

Tailored references for the individual operation modes can be accessed
via ``snekmer {mode} --help``.

.. _getting_started-configuration:

Configuration
-------------

To run Snekmer, create a ``config.yaml`` file containing desired
parameters. A `template <https://github.com/PNNL-CompBio/Snekmer/blob/main/resources/config.yaml>`_
is included in the repository. Note that a config file must be
included, in the same directory as input directory, for Snekmer
to operate.

Snekmer assumes that input files are stored in the ``input`` directory,
and automatically creates an ``output`` directory to save all output
files. Snekmer also assumes background files, if any, are stored in
``input/background``. An example of the assumed directory structure
is shown below:


Snekmer ``cluster``, ``model``, ``search``, and ``motif`` input

.. code-block:: console

    .
    ├── config.yaml
    ├── input/
    │   ├── background/
    │   │   ├── X.fasta
    │   │   ├── Y.fasta
    │   │   └── etc.
    │   ├── A.fasta
    │   ├── B.fasta
    │   └── etc.
    ├── output/
    │   ├── ...
    │   └── ...
  
  
Snekmer ``learn`` input

.. code-block:: console

    .
    ├── config.yaml
    ├── input/
    │   ├── A.fasta # known sequences to "learn" kmer counts matrix from
    │   ├── B.fasta # known sequences to "learn" kmer counts matrix from
    │   └── etc.
    │   └── base/  # optional
    │      └── base-kmer-counts.csv # optional file to additively merge kmer counts with
    ├── annotations/
    │   └── annotations.ann # annotation files used for predicting future sequences
    ├── output/
    │   ├── ...
    │   └── ...
  
  
Snekmer ``apply`` input

.. code-block:: console

    .
    ├── config.yaml
    ├── input/
    │   ├── A.fasta # unknown sequences to "apply" kmer counts matrix on
    │   ├── B.fasta # unknown sequences to "apply" kmer counts matrix on
    │   └── etc.
    ├── counts/
    │   └── kmer-counts-total.csv #kmer counts matrix generated in ``learn``
    ├── confidence/
    │   └── global-confidence-scores.csv #global confidence distribution generated in ``learn``
    ├── output/
    │   ├── ...
    │   └── ...


Partial Workflow
----------------

To execute only a part of the workflow, the ``--until`` option can be invoked.
For instance, to execute the workflow only through the kmer vector generation
step, run:

.. code-block:: bash

    snekmer {mode} --until vectorize

.. _getting_started-all_options:

All Options
-----------

.. argparse::
   :module: snekmer.cli
   :func: get_main_args
   :prog: snekmer