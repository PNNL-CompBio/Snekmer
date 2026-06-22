Command Line Interface
======================

To run any Snekmer operation mode, call:

.. code-block:: bash

    snekmer {mode}

where ``{mode}`` is one of ``easy``, ``learn``, ``apply``,
``cluster``, ``model``, or ``search``.

``easy`` is the recommended entry point for new users. It runs
the full Learn/Apply pipeline with a single command and no manual directory
setup. See the :doc:`Quick Start <quickstart>` to get started immediately.

General usage follows the pattern:

.. code-block:: bash

    snekmer <mode> [snakemake arguments] [snekmer parameter overrides]

Snakemake arguments are passed through directly to Snakemake (they are not
Snekmer-specific). Snekmer parameters can be provided via ``config.yaml`` /
``--configfile``, or overridden via Snekmer parameter flags on the command line.

For an overview of Snekmer usage, reference the help command (``snekmer --help``).

.. code-block:: console

    $ snekmer --help
    Snekmer 1.4.1 — protein sequence fingerprinting via amino acid reduction.

    Usage:
      snekmer <mode> [options]
      snekmer <mode> --help

    Modes:
      cluster            Unsupervised clustering workflow.
      model              Train supervised models + cross-validation reports.
      search             Score sequences against trained models.
      learn              Build annotation-associated k-mer distributions + confidence evaluation.
      apply              Predict annotations using outputs from learn.
      easy               Guided front-end that runs learn then apply end-to-end.

    Global options (accepted by all modes):
      --k N           K-mer length (default: 8).
      --alphabet N    Reduced alphabet encoding 0-5 or name (default: 2 = solvacc).
      --cores N       CPU cores to use (default: all).
      --dry-run       Show what would be done without executing.
      --configfile    Path(s) to YAML/JSON config file(s).
      -v, --version   Print version and exit.
      -h, --help      Show this help message and exit.

    Run 'snekmer <mode> --help' for full options for a specific mode.

Tailored references for the individual operation modes can be accessed
via ``snekmer {mode} --help``. Each subcommand includes only the Snekmer
parameter sections relevant to that mode.


.. _getting_started-configuration:

.. raw:: html

   <details class="snekmer-collapsible">
   <summary>Configuration</summary>
   <div class="details-body">

Config Precedence
`````````````````

Snekmer resolves configuration using the following precedence order (lowest to
highest):

1. **Default configfile (auto):** ``./config.yaml`` (or ``<DIR>/config.yaml``
   when using ``-d``/``--directory``).
2. **Explicit configfiles:** Any ``--configfile PATH`` values, applied in the
   order given.
3. **Snekmer parameter flags:** Any Snekmer-specific flags you explicitly
   provide on the command line (e.g. ``--k 10``, ``--alphabet hydro``).
4. **Key=Value overrides:** Any ``-C``/``--config KEY=VALUE`` overrides
   (highest precedence).

The defaults shown for Snekmer parameter flags match the template
``config.yaml`` defaults. These defaults are applied automatically only when
no config file is in use, or when a flag is explicitly provided on the command
line.

Config File
```````````

To run Snekmer with a config file, create a ``config.yaml`` file containing
desired parameters. A
`template <https://github.com/PNNL-CompBio/Snekmer/blob/main/resources/config.yaml>`_
is included in the repository.

By default, Snekmer auto-loads ``./config.yaml`` (or ``<DIR>/config.yaml``
when using ``-d``/``--directory``). You can specify one or more explicit config
files with ``--configfile``, or suppress the default auto-load with
``--no-default-configfile``.

Running Without a Config File
`````````````````````````````

A config file is no longer strictly required. You can use
``--no-default-configfile`` and rely on built-in defaults, providing any
needed overrides via Snekmer parameter flags or ``-C KEY=VALUE``.


**Directory Structure**

Snekmer assumes that input files are stored in the ``input`` directory
(configurable via ``--input-dir``), and automatically creates an ``output``
directory to save all output files. Snekmer also assumes background files,
if any, are stored in ``input/background``. An example of the assumed
directory structure is shown for each execution mode of Snekmer.

**Snekmer cluster, model, and search**

.. code-block:: console

    .
    ├── config.yaml          (optional with --no-default-configfile)
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


**Snekmer learn**

.. code-block:: console

    .
    ├── config.yaml          (optional with --no-default-configfile)
    ├── input/
    │   ├── A.fasta
    │   ├── B.fasta
    │   └── etc.
    │   └── base/            (optional)
    │      └── base-kmer-counts.csv
    ├── annotations/
    │   └── annotations.ann
    ├── output/
    │   ├── ...
    │   └── ...


**Snekmer apply**

.. code-block:: console

    .
    ├── config.yaml          (optional with --no-default-configfile)
    ├── input/
    │   ├── A.fasta
    │   ├── B.fasta
    │   └── etc.
    ├── counts/
    │   └── kmer_counts_total.csv
    ├── confidence/
    │   └── global_confidence_scores.csv
    ├── stats/
    │   └── family_summary_stats.csv
    ├── output/
    │   ├── ...
    │   └── ...

.. raw:: html

   </div>
   </details>


.. _getting_started-alphabets:

.. raw:: html

   <details class="snekmer-collapsible">
   <summary>Alphabets</summary>
   <div class="details-body">

Snekmer supports several reduced amino acid alphabets for k-mer recoding.
You may pass either an integer (``0``–``5``), the alphabet name (e.g.
``hydro``), or ``None`` to the ``--alphabet`` flag.

.. list-table::
   :header-rows: 1
   :widths: 10 20 10 60

   * - ID
     - Name
     - Size
     - Description
   * - 0
     - hydro
     - 2
     - 2-value hydrophobicity alphabet
   * - 1
     - standard
     - 7
     - "Standard" reduction alphabet
   * - 2
     - solvacc
     - 3
     - Solvent accessibility alphabet
   * - 3
     - hydrocharge
     - 3
     - 2-value hydrophobicity with charged residues as a third category
   * - 4
     - hydrostruct
     - 3
     - 2-value hydrophobicity with structural-breakers as a third category
   * - 5
     - miqs
     - 10
     - MIQS alphabet
   * - None
     - None
     - 20
     - No reduced alphabet

.. raw:: html

   </div>
   </details>


.. _getting_started-example-learnapp:

.. raw:: html

   <details class="snekmer-collapsible">
   <summary>Example: Learn→Apply Without a Config File</summary>
   <div class="details-body">

The following walkthrough demonstrates a complete ``learn`` then ``apply``
workflow using only command line arguments with no ``config.yaml`` required.
The ``--no-default-configfile`` flag tells Snekmer to skip auto-loading a
config file, so all parameters come from built-in defaults and any explicit
CLI flags.

**Step 1: Prepare the learn directory**

Create a working directory for the learn step with the expected layout:

.. code-block:: bash

    mkdir -p learn/input learn/annotations

Copy your training FASTA files and annotation file into place:

.. code-block:: bash

    cp training_sequences_*.fasta learn/input/
    cp annotations.ann            learn/annotations/

Your directory should look like:

.. code-block:: console

    learn/
    ├── annotations/
    │   └── annotations.ann
    └── input/
        ├── training_sequences_1.fasta
        ├── training_sequences_2.fasta
        └── ...

**Step 2: Run snekmer learn**

.. code-block:: bash

    snekmer learn \
        --no-default-configfile \
        --k 8 \
        --alphabet 2 \
        --input-dir input \
        --input-file-exts fasta fna faa fa \
        --input-file-regex ".*" \
        --no-nested-output \
        --no-save-apply-associations \
        --conf-weight-modifier 20 \
        --selection top_hit \
        --threshold Median \
        --apply-output snekmer_results.csv \
        -d learn

.. note::

   The values shown above match the built-in defaults. In practice, you only
   need to pass ``--no-default-configfile`` plus whichever parameters you want
   to change. For example, a minimal invocation relying entirely on defaults:

   .. code-block:: bash

       snekmer learn --no-default-configfile -d learn

   Advanced options such as sequence fragmentation are available via
   ``config.yaml``. See :doc:`config` for the full parameter reference.

**Step 3: Copy learn outputs into the apply directory**

After ``learn`` completes, create the ``apply`` directory and copy the
handoff files:

.. code-block:: bash

    mkdir -p apply/input apply/counts apply/confidence apply/stats

    cp test_sequences.fasta                                      apply/input/

    cp learn/apply_inputs/counts/kmer_counts_total.csv           apply/counts/
    cp learn/apply_inputs/confidence/global_confidence_scores.csv apply/confidence/
    cp learn/apply_inputs/stats/family_summary_stats.csv         apply/stats/

Your ``apply`` directory should look like:

.. code-block:: console

    apply/
    ├── confidence/
    │   └── global_confidence_scores.csv
    ├── counts/
    │   └── kmer_counts_total.csv
    ├── input/
    │   └── test_sequences.fasta
    └── stats/
        └── family_summary_stats.csv

**Step 4: Run snekmer apply**

.. code-block:: bash

    snekmer apply \
        --no-default-configfile \
        --k 8 \
        --alphabet 2 \
        --input-dir input \
        --input-file-exts fasta fna faa fa \
        --input-file-regex ".*" \
        --no-nested-output \
        --no-save-apply-associations \
        --conf-weight-modifier 20 \
        --selection top_hit \
        --threshold Median \
        --apply-output snekmer_results.csv \
        -d apply

.. important::

   Use the **same** ``--k`` and ``--alphabet`` values for both ``learn`` and
   ``apply``. Mismatched encoding parameters will produce incorrect results.

**Step 5: Inspect results**

The final predictions are written to ``apply/snekmer_results.csv``. You can
preview them with:

.. code-block:: bash

    head apply/snekmer_results.csv

.. raw:: html

   </div>
   </details>


.. _getting_started-partial-workflow:

.. raw:: html

   <details class="snekmer-collapsible">
   <summary>Partial Workflow</summary>
   <div class="details-body">

To execute only a part of the workflow, the ``--until`` option can be invoked.
For instance, to execute the workflow only through the kmer vector generation
step, run:

.. code-block:: bash

    snekmer {mode} --until vectorize

.. raw:: html

   </div>
   </details>


.. _getting_started-snakemake-passthrough:

.. raw:: html

   <details class="snekmer-collapsible">
   <summary>Snakemake Pass-Through Arguments</summary>
   <div class="details-body">

The following arguments are passed through directly to Snakemake and are
not Snekmer-specific:

``-n``, ``--dry-run``, ``--dryrun``
    Do not execute anything; display what would be done.

``--configfile PATH [PATH ...]``
    Specify or overwrite workflow config file(s). Multiple files overwrite
    each other in the given order.

``-C``, ``--config KEY=VALUE [KEY=VALUE ...]``
    Set or overwrite values in the workflow config object.

``--unlock``
    Unlock the working directory.

``-U``, ``--until TARGET [TARGET ...]``
    Run the workflow until the specified rules or files.

``-k``, ``--keepgoing``, ``--keep-going``
    Continue with independent jobs if a job fails.

``-w``, ``--latency``, ``--latency-wait``, ``--output-wait`` SECONDS
    Wait given seconds for output files to appear after job completion
    (default: 30).

``-t``, ``--touch``
    Touch output files instead of running commands.

``-c``, ``--cores`` N
    Use at most N CPU cores/jobs in parallel (default: all available).

``--count`` N
    Number of files to process (limits DAG size).

``--countstart`` IDX
    Starting file index for use with ``--count`` (default: 0).

``--verbose``
    Show additional debug output.

``-q``, ``--quiet`` [progress|rules|all]
    Reduce Snakemake output.

``-d``, ``--directory`` DIR
    Specify working directory.

``-R``, ``--forcerun`` [TARGET ...]
    Force re-execution/creation of the given rules or files.

``--list-code-changes``, ``--lc``
    List output files for which the rule body changed.

``--list-params-changes``, ``--lp``
    List output files for which defined params changed.

``--no-default-configfile``
    Do not auto-load ``./config.yaml`` (or ``<DIR>/config.yaml`` with ``-d``).

``--clust`` PATH [PATH ...]
    Path to cluster execution YAML configuration file (e.g., for SLURM).

``-j``, ``--jobs`` N
    Number of simultaneous jobs to submit to the scheduler (default: 1000).

``--scheduler`` **[greedy|ilp]**
    Specify whether Snakemake uses the greedy or ILP scheduler.

``--scheduler-ilp-solver`` **[SOLVER]**
    Specify the MILP solver to be used when using the ILP scheduler.

``--scheduler-ilp-solver-path`` **[PATH]**
    PATH to search for ILP scheduler solver binaries.

.. raw:: html

   </div>
   </details>


.. _getting_started-all_options:

.. raw:: html

   <details class="snekmer-collapsible">
   <summary>All Options (full argparse reference)</summary>
   <div class="details-body">

.. argparse::
   :module: snekmer.cli
   :func: get_main_args
   :prog: snekmer

.. raw:: html

   </div>
   </details>
