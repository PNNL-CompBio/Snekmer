Command Line Interface
======================

To run any of the 3 Snekmer operation modes, simply call:

.. code-block:: bash

    snekmer {mode}

Each mode has its own mode-specific options and parameters to be specified
on the command line or the ``config.yaml`` file, respectively.

For an overview of Snekmer usage, reference the help command (``snekmer --help``).

.. code-block:: console

    $ snekmer --help
    usage: snekmer [-h] [-v] [--dryrun] [--configfile PATH] [--config [KEY=VALUE ...]] [--unlock] [--until TARGET [TARGET ...]]
               [--latency SECONDS] [--touch] [--cores N] [--count N] [--countstart IDX] [--cluster PATH] [--jobs N]
               {cluster,model,search}

    Snekmer: A tool for kmer-based sequence analysis using amino acid reduction (AAR)

    positional arguments:
    {cluster,model,search}

    options:
    -h, --help            show this help message and exit
    -v, --version         print version and exit
    --dryrun              perform a dry run
    --configfile PATH     path to yaml configuration file
    --config [KEY=VALUE ...]
                            Set or overwrite values in the workflow config object. The workflow config object is accessible as
                            variable config inside the workflow. Default values can be set by providing a JSON file.
    --unlock              unlock directory
    --until TARGET [TARGET ...]
                            run pipeline until reaching the target rule or files
    --latency SECONDS, --latency-wait SECONDS
                            wait time, in seconds, for output file creation (default 30)
    --touch               touch output files only
    --cores N             number of cores used for execution (local execution only)
    --count N             number of files to process (limits DAG size)
    --countstart IDX      starting file index (for use with --count)

    cluster arguments:
    --cluster PATH        path to cluster execution yaml configuration file
    --jobs N              number of simultaneous jobs to submit to a slurm queue

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

Partial Workflow
----------------

To execute only a part of the workflow, the ``--until`` option can be invoked.
For instance, to execute the workflow only through the kmer vector generation
step, run:

.. code-block:: bash
    
    snekmer {mode} --until vectorize
