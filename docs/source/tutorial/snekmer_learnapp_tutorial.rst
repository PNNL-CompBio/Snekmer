.. _learnapp-tutorial:

Learn/Apply — Full Pipeline Reference
======================================

.. note::

   **New users:** ``snekmer easy-learn-apply`` is the recommended entry point.
   It runs the complete pipeline with a single command and no manual setup.
   See the :doc:`easy-learn-apply tutorial <snekmer_easy_learn_apply_tutorial>`.

   This page is for users who need **direct control** over ``snekmer learn``
   and ``snekmer apply`` — for example:

   - Incrementally adding new training sequences to an existing model
   - Reusing a trained model against multiple query sets
   - Customising intermediate pipeline steps


When to use ``learn`` / ``apply`` directly
-------------------------------------------

.. list-table::
   :header-rows: 1
   :widths: 55 45

   * - Situation
     - Recommendation
   * - First time running, want results fast
     - ``easy-learn-apply``
   * - Existing ``config.yaml`` and directory layout
     - ``snekmer learn`` then ``snekmer apply``
   * - Adding new training data to an existing model
     - ``snekmer learn`` then ``snekmer apply``
   * - Reusing a trained model against new query sequences
     - ``snekmer apply`` only (skip learn)


Demo data
---------

The commands below use the demo data in ``resources/demo_sequences/learn_apply_inputs/``:

.. code-block:: text

   resources/demo_sequences/learn_apply_inputs/
   ├── learn/                        ← 10 training FASTA files (10,000 proteins, 200 families)
   │   ├── training_sequences_1.fasta
   │   ├── ...
   │   └── training_sequences_10.fasta
   ├── apply/
   │   └── test_sequences_1.fasta    ← 3,000 query proteins
   └── annotations/
       └── TIGRFAMs_annotation.ann   ← id/family TSV

All commands assume you are running from the **root of the Snekmer repository**.


Directory layout
----------------

``snekmer learn`` and ``snekmer apply`` each require their own working directory
with a specific structure. ``easy-learn-apply`` creates these automatically;
when using the modes directly you build them yourself.

``learn`` workspace
~~~~~~~~~~~~~~~~~~~

.. code-block:: text

   learn/
   ├── config.yaml
   ├── annotations/
   │   └── annotations.ann       ← tab-separated id / family file
   └── input/
       ├── training_sequences_1.fasta
       └── ...

``apply`` workspace
~~~~~~~~~~~~~~~~~~~

.. code-block:: text

   apply/
   ├── config.yaml
   ├── input/
   │   └── test_sequences_1.fasta
   ├── counts/
   │   └── kmer_counts_total.csv      ← copied from learn output
   ├── confidence/
   │   └── global_confidence_scores.csv   ← copied from learn output
   └── stats/
       └── family_summary_stats.csv   ← copied from learn output


Step 1 — Set up the ``learn`` workspace
----------------------------------------

.. code-block:: bash

   mkdir -p learn/input learn/annotations

   # Copy training sequences into the workspace
   cp resources/demo_sequences/learn_apply_inputs/learn/*.fasta learn/input/

   # Copy annotation file
   cp resources/demo_sequences/learn_apply_inputs/annotations/TIGRFAMs_annotation.ann \
      learn/annotations/annotations.ann

.. note::

   On Linux/macOS you can use symlinks instead of copying to save disk space:

   .. code-block:: bash

      ln -s "$(pwd)/resources/demo_sequences/learn_apply_inputs/learn/"*.fasta learn/input/

   Symlinks are not supported on Windows — use ``cp`` there.


Step 2 — Run ``snekmer learn``
--------------------------------

If you have a ``config.yaml`` in the ``learn/`` directory:

.. code-block:: bash

   snekmer learn -d learn

Without a config file (uses built-in defaults):

.. code-block:: bash

   snekmer learn --no-default-configfile -d learn

To preview what will run without executing:

.. code-block:: bash

   snekmer learn --no-default-configfile --dry-run -d learn

.. tip::

   Use the same ``--k`` and ``--alphabet`` values for both ``learn`` and ``apply``.
   Mismatched encoding parameters will produce incorrect results.

``learn`` writes its outputs to ``learn/output/`` and creates a convenience
``learn/apply_inputs/`` directory alongside it:

.. code-block:: text

   learn/
   ├── apply_inputs/           ← ready-to-use handoff files for snekmer apply
   │   ├── counts/kmer_counts_total.csv
   │   ├── stats/family_summary_stats.csv
   │   └── confidence/global_confidence_scores.csv
   └── output/
       ├── kmerize/    ← per-file k-mer labels (.kmers)
       ├── vector/     ← per-file k-mer vectors (.npz)
       ├── learn/      ← per-file and merged k-mer count matrices
       └── eval_conf/  ← confidence scores and family statistics


Step 3 — Copy ``learn`` outputs into the ``apply`` workspace
--------------------------------------------------------------

.. code-block:: bash

   mkdir -p apply/input apply/counts apply/confidence apply/stats

   # Query sequences (use cp on Windows; symlink on Linux/macOS)
   cp resources/demo_sequences/learn_apply_inputs/apply/test_sequences_1.fasta apply/input/

   # Handoff files from learn (apply_inputs/ is at the root of the learn workspace)
   cp learn/apply_inputs/counts/kmer_counts_total.csv       apply/counts/
   cp learn/apply_inputs/confidence/global_confidence_scores.csv  apply/confidence/
   cp learn/apply_inputs/stats/family_summary_stats.csv     apply/stats/


Step 4 — Run ``snekmer apply``
--------------------------------

.. code-block:: bash

   snekmer apply -d apply

Or without a config file:

.. code-block:: bash

   snekmer apply --no-default-configfile -d apply

Results are written to ``apply/snekmer_results.csv`` (one row per query sequence):

.. code-block:: text

   Sequence                          Prediction   Score    delta   Confidence
   tr|A0A2S8EUS7|A0A2S8EUS7_9RHOB   TIGR01783    0.199    0.00    0.383
   tr|A0A401ZGP4|A0A401ZGP4_9CHLR   TIGR00757    0.316    0.02    0.922
   tr|A0A427BXE3|A0A427BXE3_9GAMM   TIGR01023    0.198    0.08    1.000
   ...

See :doc:`snekmer_easy_learn_apply_tutorial` for a description of each output column
and guidance on filtering by confidence score.


Key parameters
--------------

The most commonly adjusted parameters. Pass them as CLI flags or set them in
``config.yaml`` — see :doc:`../getting_started/config` for the full reference.

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   * - Parameter
     - Default
     - Description
   * - ``--k`` / ``k``
     - ``8``
     - K-mer length
   * - ``--alphabet`` / ``alphabet``
     - ``2`` (solvacc)
     - Amino acid reduction alphabet (0–5 or name)
   * - ``--selection`` / ``selection``
     - ``top_hit``
     - Annotation selection method: ``top_hit``, ``greatest_distance``, ``combined_distance``
   * - ``--threshold`` / ``threshold``
     - ``Median``
     - Family score threshold for filtering: ``Median``, ``Mean``, ``90th Percentile``, ``None``
   * - ``--apply-output`` / ``apply_output``
     - ``snekmer_results.csv``
     - Output filename for the results CSV

For the full list of options run ``snekmer learn --help`` or ``snekmer apply --help``,
or see :ref:`All Options <getting_started-all_options>` in the CLI reference.


Reusing a trained model
-----------------------

If you already have a ``learn/`` workspace with valid ``apply_inputs/`` (e.g., from a
previous run), you can skip ``snekmer learn`` and run ``snekmer apply`` against any
new set of query sequences. Just update ``apply/input/`` and re-run apply.

To add new training families to an existing model, re-run ``snekmer learn`` pointing
to the expanded training set. The merged counts matrix accumulates across runs.


Deep-dive notebook
------------------

For a step-by-step walkthrough of every internal pipeline rule — vectorization,
k-mer count matrix construction, reversed-sequence decoy evaluation, and confidence
calibration — see the companion notebook:

``docs/source/tutorial/snekmer_learn_apply_tutorial.ipynb``

This notebook exposes the Python code behind each Snakemake rule and is intended
for users who want to understand the method in detail or adapt intermediate outputs.
