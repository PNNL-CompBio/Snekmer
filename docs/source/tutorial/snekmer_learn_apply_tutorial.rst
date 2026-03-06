========================================
Snekmer Learn/Apply Tutorial
========================================

This tutorial walks through the process of running the **Snekmer Learn** and
**Snekmer Apply** pipelines, from directory setup through vectorization and
annotation prediction.

Overview
========

The basic process for running Snekmer Learn/Apply is as follows:

1. Verify that your file directory structure is correct and that the top-level
   directory contains a ``config.yaml`` file.

   - A ``config.yaml`` template is included in the Snekmer codebase at
     ``resources/learn_apply/config.yaml``.

2. Modify ``config.yaml`` with the desired parameters.

3. Use the command line to navigate to the directory containing both the
   ``config.yaml`` file and the ``input/`` directory.

4. Run ``snekmer learn``, then copy the appropriate outputs to a separate
   directory to run ``snekmer apply``.


Running the Snekmer Learn Pipeline
===================================

Setup
-----

.. note::

   This tutorial assumes you are working from the ``resources/tutorial/``
   directory.

To set up the workflow we initialize a configuration dictionary (equivalent to
the YAML file used on the command line) and gather all input files.  Input files
are detected with ``glob.glob``, exactly as Snekmer performs input-file
detection.

.. code-block:: python

   # --- Standard library ---
   import os, sys, shutil, gzip, pickle
   from glob import glob
   from pathlib import Path

   # --- Third-party ---
   # (additional imports as required by your environment)


Configuration Parameters
~~~~~~~~~~~~~~~~~~~~~~~~~

Below is an example configuration dictionary.  When using the command-line
workflow these values live in ``config.yaml``.

.. code-block:: python

   config = {
       # ...
       "save_apply_associations": False,

       # Fragmentation settings
       "fragmentation": False,
       "version": "absolute",
       "frag_length": 50,
       "min_length": 50,
       "location": "random",
       "seed": 999,

       # Confidence weighting
       "conf_weight_modifier": 20,

       # Selection methods: "top_hit", "greatest_distance", "combined_distance"
       "selection": "top_hit",

       # Threshold column from family_summary_stats.csv (or None)
       "threshold": "Median",

       # Weights (only relevant when using combined_distance)
       "weight_top": 0.5,
       "weight_distance": 0.5,

       # Output naming
       "apply_output": "snekmer_results.csv",
   }

.. warning::

   The ``fragmentation`` option and the ``conf_weight_modifier`` parameter
   (intended for subsequent data additions) are **not** supported in the
   notebook workflow.


Rule 0 — Collect Input Files
-----------------------------

Gather all FASTA files from the input directory and define the output
directory:

.. code-block:: python

   input_dir = "learn/input"

   unzipped = []
   for ext in ["fasta", "fa", "faa"]:
       unzipped.extend(glob(os.path.join(input_dir, f"*.{ext}")))
   print("unzipped files:\t", unzipped)

   # Define output directory (create if missing)
   output_dir = "learn/output"
   os.makedirs(output_dir, exist_ok=True)


Vectorize Helper Function
--------------------------

The helper function below replicates the k-mer vectorization step normally
performed inside the Snakemake workflow (``kmerize``, ``vectorize``).  Given a
FASTA file it:

* constructs (or loads) a k-mer basis,
* reduces amino-acid sequences to the chosen alphabet,
* generates a binary presence/absence matrix of k-mers for every protein, and
* writes the results to the ``.npz`` and ``.kmers`` formats used by the full
  Learn pipeline.

.. code-block:: python

   def run_vectorize_like_snakemake(
       fasta_path: str,
       alphabet,
       k: int,
       output_npz: str,
       output_kmerobj: str,
       basis_path: str = None,
       min_filter: int = 0,
   ):
       """Notebook helper that replicates kmerize.smk::vectorize."""
       kmer = KmerVec(alphabet=alphabet, k=k)
       if basis_path is not None and os.path.exists(basis_path):
           # ... load existing basis and vectorize ...
           pass
       # (full implementation omitted for brevity)


Vectorize Input Files
---------------------

Create the output directories and vectorize every FASTA file collected in
**Rule 0**:

.. code-block:: python

   vector_dir = os.path.join(output_dir, "vector")
   kmer_dir   = os.path.join(output_dir, "kmerize")
   os.makedirs(vector_dir, exist_ok=True)
   os.makedirs(kmer_dir,   exist_ok=True)

   for fa in unzipped:
       base    = os.path.basename(skm.utils.split_file_ext(fa)[0])
       out_npz = os.path.join(vector_dir, f"{base}.npz")
       # ... call run_vectorize_like_snakemake(...) ...

Expected output:

.. code-block:: text

   Vectorized learn/input/training_sequences_1.fasta  → learn/output/vector/training_sequences_1.npz
   Vectorized learn/input/training_sequences_2.fasta  → learn/output/vector/training_sequences_2.npz
   ...
   Vectorized learn/input/training_sequences_10.fasta → learn/output/vector/training_sequences_10.npz


Visualising Prediction Results
==============================

After running both **Learn** and **Apply**, the annotation predictions can be
visualised with a simple bar chart:

.. code-block:: python

   fig, ax = plt.subplots()
   # ... (build bar chart from results DataFrame) ...
   ax.set_title("Annotation Prediction Results")
   ax.set_ylabel("Sequences")
   h, l = ax.get_legend_handles_labels()
   u = dict(zip(l, h))
   ax.legend(u.values(), u.keys(), ncol=2, fontsize=8)
   plt.tight_layout()
   plt.show()
