Using Snekmer
=============

Annotation Pipeline (Learn/Apply)
-----------------------------------

The primary use case for Snekmer is sequence annotation via the Learn/Apply pipeline,
accessible through three commands:

- ``easy-learn-apply``: **recommended entry point** that runs the full Learn/Apply pipeline
  from a single command. Provide training sequences, query sequences, and an annotation file;
  Snekmer handles the rest:

  .. code-block:: bash

     snekmer easy-learn-apply --train train/ --query query.fasta --ann annotations.ann --output results/

- ``learn``: builds a kmer counts matrix and confidence model from annotated training sequences.
- ``apply``: scores query sequences against outputs from a prior ``learn`` run.

The ``easy-learn-apply`` command is built on top of ``learn`` and ``apply`` and produces
identical results. Use ``learn``/``apply`` directly when adding new training data to an existing
model or when you need fine-grained control over intermediate steps.


Additional Modes (Cluster / Model / Search)
--------------------------------------------

Snekmer also supports unsupervised clustering, supervised modeling, and model-based search.
See the :doc:`Model/Cluster/Search tutorial <../tutorial/snekmer_demo>` for a full walkthrough.

These modes share a common directory layout:

.. code-block:: text

   my_project/
   ├── config.yaml      ← copy from resources/config.yaml and edit
   └── input/
       ├── family_A.fasta
       ├── family_B.fasta
       └── ...          ← one FASTA file per protein family

Each FASTA file should contain sequences belonging to a single protein family.
The filename (without extension) is used as the family label.

Run from the ``my_project/`` directory:

.. code-block:: bash

    snekmer cluster    # unsupervised clustering
    snekmer model      # supervised ML models (one-vs-rest)
    snekmer search     # score unknowns against trained models

Preview the pipeline steps without executing with ``--dryrun``:

.. code-block:: bash

    snekmer model --dryrun

An example ``config.yaml`` is included at
`resources/config.yaml <https://github.com/PNNL-CompBio/Snekmer/blob/main/resources/config.yaml>`_.

.. _usage-results:

Accessing Results
-----------------

Summary Reports
:::::::::::::::

Each step in the Snekmer modeling pipeline will generate a report in HTML format.
Users can find these reports, entitled **Snekmer_\<MODE\>_Report.html**,
in the output directory.

Common Output Files (all modes)
::::::::::::::::::::::::::::::::

All operation modes preprocess input files and kmerize sequences.
The associated output files can be found in the respective directories.

The following output directories and files will always be created:

.. code-block:: console

    .
    ├── input/
    │   ├── A.fasta
    │   └── B.fasta
    ├── output/
    │   ├── kmerize/
    │   │   ├── A.kmers  # kmer labels for A
    │   │   └── B.kmers  # kmer labels for B
    │   ├── vector/
    │   │   ├── A.npz    # sequences, sequence IDs, and kmer vectors for A
    │   │   └── B.npz    # sequences, sequence IDs, and kmer vectors for B
    │   ├── ...

Mode-Specific Output Files
--------------------------

The steps in the Snekmer pipeline generate their own associated output files.

Snekmer Cluster Output Files
::::::::::::::::::::::::::::

Snekmer's cluster mode produces the following output files
and directories in addition to the files described previously.

.. code-block:: console

    .
    └── output/
        ├── ...
        └── cluster/
            ├── snekmer.csv     # Summary of clustering results
            └── figures/        # Clustering figures
                ├── pca_explained_variance_curve.png
                ├── tsne.png
                └── umap.png

Snekmer Model Output Files
::::::::::::::::::::::::::

Snekmer's model mode produces the following output files
and directories in addition to the files described previously.

.. code-block:: console

    .
    ├── output/
    │   ├── ...
    │   ├── scoring/
    │   │   ├── A.matrix    # Similarity matrix for A seqs
    │   │   ├── B.matrix    # Similarity matrix for B seqs
    │   │   ├── A.scorer    # Object to apply A scoring model
    │   │   ├── B.scorer    # Object to apply B scoring model
    │   │   └── weights/
    │   │       ├── A.csv.gz    # Kmer score weights in A kmer space
    │   │       └── B.csv.gz    # Kmer score weights in B kmer space
    │   ├── model/
    │   │   ├── A.model     # (A/not A) classification model
    │   │   ├── B.model     # (B/not B) classification model
    │   │   ├── results/    # Cross-validation results tables
    │   │   │   ├── A.csv
    │   │   │   └── B.csv
    │   │   └── figures/      # Cross-validation results figures
    │   │       ├── A/
    │   │       └── B/

Snekmer Search Output Files
:::::::::::::::::::::::::::

The ``snekmer search`` mode assumes that the user has pre-generated
family models using the ``snekmer model`` workflow, and thus operates
as an independent workflow. The location of the basis sets, scorers,
and models must be specified in the configuration file (see the search
params section in the provided
`example <https://github.com/PNNL-CompBio/Snekmer/blob/main/resources/config.yaml>`_).

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

The user should then modify their configuration file to point towards
the appropriate basis set, scorer, and model directories.

Executing ``snekmer search --configfile search.yaml`` produces the
following output files and directories in addition to the files
described previously.

.. code-block:: console

    .
    └── output/
        ├── kmers/
        │   └── common.basis  # Common kmer basis set for queried families
        └── search/
            ├── A   # A probabilities and predictions for unknown sequences
            │   ├── unknown_1.csv
            │   ├── unknown_2.csv
            │   └── ...
            └── B   # B probabilities and predictions for unknown sequences
                ├── unknown_1.csv
                ├── unknown_2.csv
                └── ...  



Snekmer Learn Output Files
::::::::::::::::::::::::::

Snekmer's learn mode produces the following output files
and directories in addition to the files described previously.

.. code-block:: console

    .
    ├── apply_inputs/           ← ready-to-use handoff files for snekmer apply
    │   ├── counts/
    │   │   └── kmer_counts_total.csv
    │   ├── stats/
    │   │   └── family_summary_stats.csv
    │   └── confidence/
    │       └── global_confidence_scores.csv
    └── output/
        ├── kmerize/
        │   ├── A.kmers  # kmer labels for A
        │   └── B.kmers  # kmer labels for B
        ├── vector/
        │   ├── A.npz    # sequences, sequence IDs, and kmer vectors for A
        │   └── B.npz    # sequences, sequence IDs, and kmer vectors for B
        ├── learn/
        │   ├── kmer_counts_A.csv        # kmer counts matrix for A seqs
        │   ├── kmer_counts_B.csv        # kmer counts matrix for B seqs
        │   └── kmer_counts_total.csv    # merged kmer counts matrix
        ├── eval_conf/
        │   ├── global_confidence_scores.csv    # global confidence score distribution
        │   ├── family_summary_stats.csv        # per-family score statistics
        │   └── family_stats_checkpoint.csv     # incremental-update checkpoint
        └── evaluate/
            ├── eval_apply_sequences/
            │   └── seq-annotation-scores-A.model   # self-assessed cosine similarity scores
            └── eval_apply_reversed/
                └── seq-annotation-scores-A.csv.gz  # scores for reversed-sequence decoys

Snekmer Apply Output Files
::::::::::::::::::::::::::

Snekmer's apply mode produces the following output files
and directories in addition to the files described previously.
Predictions are stored in ``kmer_summary_<name>.csv``, one 5-column file per input FASTA,
with one row per sequence: the predicted family, cosine similarity score, score gap (delta),
and calibrated confidence. These are concatenated into the single ``snekmer_results.csv`` file.
The optional ``seq_annotation_scores_<name>.csv`` files contain the full cosine similarity
matrix (one row per sequence, one column per family) and can be large for big datasets.

.. code-block:: console

    .
    ├── snekmer_results.csv          ← compiled predictions (all sequences)
    └── output/
        ├── apply/
        │   ├── kmer_summary_C.csv              ← predictions and confidence for C seqs
        │   ├── kmer_summary_D.csv              ← predictions and confidence for D seqs
        │   ├── seq_annotation_scores_C.csv     ← (optional) all cosine similarity scores for C
        │   └── seq_annotation_scores_D.csv     ← (optional) all cosine similarity scores for D
        └── Snekmer_Apply_Report.html


.. _usage-interpreting:

Interpreting Predictions
------------------------

The main output of any Learn/Apply run is ``snekmer_results.csv``. It contains
one row per query sequence with five columns:

.. list-table::
   :header-rows: 1
   :widths: 18 82

   * - Column
     - Meaning
   * - ``Sequence``
     - Sequence identifier extracted from the FASTA header. For UniProt-style
       headers (``>db|ACCESSION|name``), this is the field between the first
       pair of ``|`` characters.
   * - ``Prediction``
     - The family assigned the highest cosine similarity score. Every sequence
       receives a prediction; use ``Score`` and ``Confidence`` to decide
       whether to accept it.
   * - ``Score``
     - Cosine similarity between the query k-mer vector and the predicted
       family's k-mer profile (0–1). Higher values indicate greater overlap
       with the training sequences. A score of **0.0** means the query shares
       no k-mers with any training family; these predictions are not
       meaningful and should be excluded.
   * - ``delta``
     - Gap between the top score and the second-best family score. A larger
       gap means the top assignment is more distinct from competing families.
       Low delta (near 0) indicates ambiguity between families.
   * - ``Confidence``
     - Calibrated probability that the prediction is correct (0–1), derived
       by comparing the query's score to the distribution of self-scores from
       training sequences. Values near **1.0** indicate the score is
       consistent with confidently annotated training sequences; values near
       **0.0** indicate the score is below the typical range seen in training.

**Recommended filtering**

A ``Confidence ≥ 0.95`` threshold is a good starting point for high-quality
annotations:

.. code-block:: python

   import pandas as pd
   df = pd.read_csv("snekmer_results.csv")
   high_conf = df[(df["Confidence"] >= 0.95) & (df["Score"] > 0)]
   print(f"{len(high_conf)} high-confidence predictions out of {len(df)}")
   high_conf.to_csv("high_confidence_results.csv", index=False)

Adjusting the threshold trades precision for recall:

- **Higher threshold** (e.g. 0.99): fewer predictions, greater reliability.
- **Lower threshold** (e.g. 0.80): more predictions, more false positives.

**Low-confidence predictions**

Sequences with low ``Confidence`` or ``Score = 0`` may:

- Belong to families not represented in the training set.
- Be too divergent for the chosen k-mer length or alphabet.
- Indicate that more training sequences are needed for those families
  (at least 20–50 sequences per family is recommended).

