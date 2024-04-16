Using Snekmer
=============

Snekmer has three modeling operations: ``cluster`` (unsupervised clustering),
``model`` (supervised modeling), and ``search`` (application
of model to new sequences). We will call the first two modes
**learning modes** due to their utility in learning relationships
between protein family input files. Users may choose a mode to best
suit their specific use case.

Snekmer also has two non-modeling operations: ``learn`` (kmer counts matrix generation), 
and ``apply`` (perform cosine similarity between sequences and kmer counts matrix). The Learn/Apply
pipeline can (and should) be used with with large training datasets to quickly find predict 
annotations for new sequences.  

The mode must be specified in the command line, e.g. to specify the
``model`` mode, the following should be called:

.. code-block:: bash

    snekmer model [--options]

In the `resources <https://github.com/PNNL-CompBio/Snekmer/tree/main/resources>`_,
an example configuration file is included:

  - `config.yaml <https://github.com/PNNL-CompBio/Snekmer/blob/main/resources/config.yaml>`_: Configuration file for snekmer execution.

.. code-block:: bash

    snekmer {mode} --dryrun

(For instance, in supervised mode, run ``snekmer model --dryrun``.)

The output of the dry run shows you the files that will be created by the
pipeline. If no files are generated, double-check   that your directory
structure matches the format specified above.

When you are ready to process your files, run:

.. code-block:: bash

    snekmer {mode}

.. _usage-results:

Accessing Results
-----------------

Summary Reports
:::::::::::::::

Each step in the Snekmer modeling pipeline will generate a report in HTML format.
Users can find these reports, entitled **Snekmer_\<MODE\>_Report.html**,
in the output directory.

Snekmer Model Output Files
::::::::::::::::::::

All operation modes will preprocess input files and kmerize sequences.
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
    ├── output/
    │   ├── kmerize/
    │   │   ├── A.kmers  # kmer labels for A
    │   │   └── B.kmers  # kmer labels for B
    │   ├── vector/
    │   │   ├── A.npz    # sequences, sequence IDs, and kmer vectors for A
    │   │   └── B.npz    # sequences, sequence IDs, and kmer vectors for B
    │   ├── vector_frag/ 
    │   │   ├── A.npz    # Conditional output for vector when the fragmentation option is True.
    │   │   └── B.npz    # Conditional output for vector when the fragmentation option is True.
    │   ├── learn/
    │   │   ├── kmer-counts-A.csv    # Kmer Counts matrix for A seqs
    │   │   ├── kmer-counts-B.csv     # Kmer Counts matrix for B seqs
    │   │   └── kmer-counts-total.csv    # Kmer Counts matrix for merged (total) database.
    │   ├── eval_apply/
    │   │   ├── Seq-Annotation-Scores-A.model     # Self-assessed sequence-annotation cosine similarity scores for A seqs
    │   │   ├── Seq-Annotation-Scores-B.model     # Self-assessed sequence-annotation cosine similarity scores for B seqs
    │   ├── eval_apply_frag/
    │   │   ├── Seq-Annotation-Scores-A.model     # Conditional output for eval_apply when the fragmentation option is True.
    │   │   ├── Seq-Annotation-Scores-B.model     # Conditional output for eval_apply when the fragmentation option is True.
    │   ├── eval_conf/
    │   │   ├── global-confidence-scores.csv     # Global confidence score distribution
    │   │   └── confidence_matrix.csv   # Confidence distribution Matrix for each annotation


Snekmer Apply Output Files
::::::::::::::::::::::::::

Snekmer's apply mode produces the following output files
and directories in addition to the files described previously.

.. code-block:: console

    .
    ├── output/
    │   ├── ...
    │   ├── apply/
    │   │   ├── Seq-Annotation-Scores-C.csv  # (optional) Sequence-annotation cosine similarity scores for C seqs
    │   │   ├── Seq-Annotation-Scores-D.csv  # (optional) Sequence-annotation cosine similarity scores for D seqs
    │   │   ├── kmer-summary-C.csv  # Results with annotation predictions and confidence for C seqs 
    │   │   └── kmer-summary-D.csv  # Results with annotation predictions and confidence for D seqs 
