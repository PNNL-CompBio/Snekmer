.. _config-main:

Setting up User Configuration (config.yaml)
===========================================

To run Snekmer, the user must specify parameters either in a configuration
file (.YAML) or by specifying them as command line arguments using -C or --config. A template ``config.yaml`` file is included in the
`resources directory <https://github.com/PNNL-CompBio/Snekmer/tree/main/resources>`_.

The example YAML files included are:

* ``config.yaml``: Configuration file for running Snekmer
* ``clust.yaml``: (optional) Cluster configuration file for deploying Snekmer on a high-performance computing (HPC) cluster

Parameter Descriptions for ``config.yaml``
------------------------------------------

The base `config.yaml` file contains the parameters which are required in order to run `snekmer model` or `snekmer cluster`. These may alternatively be specified using command line arguments, and Snekmer supports specifying some parameters in a .yaml file (either config.yaml or specified using --configfile when invoking Snekmer) and others using -C or --config arguments.

Required Parameters
```````````````````

Parameters which are required to be specified by the user in order to use Snekmer.

====================  ====================  ===================================================================================================
     Parameter                Type           Description
====================  ====================  ===================================================================================================
 ``k``                 ``int``               K-mer length
 ``alphabet``          ``str`` or ``int``    Reduced alphabet encoding
                                             (see :ref:`Alphabets <alphabets>` for more details). Alphabets may be specified by numbers 0-5 or by their names.
====================  ====================  ===================================================================================================

Input/Output Parameters
```````````````````````

General parameters related to input and output sequences and/or files.

========================  ====================  =========================================================================
     Parameter                    Type            Description
========================  ====================  =========================================================================
 ``input_dir``             ``str``               Directory containing input FASTA files (default: ``input``)
 ``input_file_exts``       ``list``               File extensions to be considered as valid for input sequence files
 ``input_file_regex``      ``str`` or ``None``    Regular expression for parsing family/annotation identifiers from filenames
 ``nested_output``         ``bool``               If True, saves files into nested directory structure, i.e. ``{save_dir}/{alphabet}/{k}``
========================  ====================  =========================================================================

Score Parameters
````````````````

General parameters related to how Snekmer calculates family scores for k-mers.

========================  =====================  =================================================================================
     Parameter                   Type             Description
========================  =====================  =================================================================================
 ``scaler``                ``bool``               If True, applies k-mer frequency scaling before scoring
 ``scaler_kwargs``         ``dict``               Keyword arguments to pass to k-mer scaler object (e.g. ``{"n": 0.25}``)
 ``labels``                ``str`` or ``None``    If None, uses default kmer set for scaler. Otherwise, uses the ones specified
 ``lname``                 ``str`` or ``None``    Label name (e.g. ``"family"``)
========================  =====================  =================================================================================

Model Parameters
````````````````

General parameters related to Snekmer's model mode (``snekmer model``), wherein supervised models are trained via the workflow.

========================  =====================  =========================================================================
     Parameter                    Type            Description
========================  =====================  =========================================================================
 ``cv``                    ``int``                Number of cross-validation folds for model evaluation
 ``random_state``          ``int`` or ``None``    Random state for model evaluation
========================  =====================  =========================================================================

Cluster Parameters
``````````````````

General parameters related to Snekmer's cluster mode (``snekmer cluster``), wherein unsupervised clusters are produced via the workflow.

========================  ====================  ==============================================================================
     Parameter                    Type            Description
========================  ====================  ==============================================================================
 ``method``                ``str``                Clustering algorithm. See table below for all options.
 ``params``                ``dict``               Parameters to pass to the clustering algorithm
 ``cluster_plots``         ``bool``               If True, generates figures illustrating clustering results (t-SNE, UMAP, PCA)
 ``min_rep``               ``int`` or ``None``    Discard k-mers with fewer than this many occurrences across the input set
 ``max_rep``               ``int`` or ``None``    Discard k-mers with more than this many occurrences across the input set
 ``save_matrix``           ``bool``               If True, saves the pairwise distance matrix (large files — not recommended for large datasets)
 ``dist_thresh``           ``int``                Distance threshold used when computing the BSF Jaccard matrix
========================  ====================  ==============================================================================

**Clustering methods** (``cluster.method``)

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Value
     - Description
   * - ``agglomerative-jaccard``
     - Agglomerative clustering using Jaccard distance (default; requires BSF or falls back to scipy)
   * - ``density-jaccard``
     - DBSCAN density clustering using Jaccard distance (requires BSF or falls back to scipy)
   * - ``hdensity-jaccard``
     - HDBSCAN density clustering using Jaccard distance (requires BSF or falls back to scipy)
   * - ``agglomerative``
     - Agglomerative clustering using Euclidean distance
   * - ``kmeans``
     - Mini-batch k-means clustering
   * - ``correlation``
     - Hierarchical clustering using correlation distance
   * - ``density``
     - DBSCAN density-based clustering
   * - ``birch``
     - Birch incremental clustering
   * - ``optics``
     - OPTICS density-based clustering
   * - ``hdbscan``
     - HDBSCAN hierarchical density clustering

The three ``-jaccard`` methods use the
:doc:`Blazing Signature Filter (BSF) <advanced>` when installed, and fall back
to ``scipy.spatial.distance.pdist`` automatically when BSF is not available.


Parameter Descriptions for ``clust.yaml``
-------------------------------------------

``clust.yaml`` is an **optional** configuration file used to deploy Snekmer jobs
on a high-performance computing (HPC) cluster via SLURM (or another scheduler
supported by Snakemake). It is not required for local runs.

A typical ``clust.yaml`` specifies resource requests per rule, for example:

.. code-block:: yaml

   __default__:
     partition: normal
     time: "04:00:00"
     mem: "16G"
     ntasks: 1
     cpus-per-task: 4

   vectorize:
     time: "01:00:00"
     mem: "8G"

Pass it to Snekmer with the ``--clust`` flag:

.. code-block:: bash

   snekmer cluster --clust clust.yaml

See the `Snakemake cluster execution documentation <https://snakemake.readthedocs.io/en/stable/executing/cluster.html>`_
and `SLURM sbatch documentation <https://slurm.schedmd.com/sbatch.html>`_ for
the full list of supported fields.

Required Parameters for Snekmer Search
--------------------------------------

The following parameters must be specified when running `snekmer search`.

========================  =====================  ========================================================================================
     Parameter                     Type           Description
========================  =====================  ========================================================================================
 ``input_file_exts``       ``list``               See `Input/Output Parameters`_
 ``input_file_regex``      ``str`` or ``None``    See `Input/Output Parameters`_
 ``model_dir``             ``str``                Directory containing model object(s) (.model)
 ``basis_dir``             ``str``                Directory containing k-mer basis set(s) (.kmers)
 ``score_dir``             ``str``                Directory containing scoring object(s) (.scorer)
 ``k``                     ``int``                See `Required Parameters`_
 ``alphabet``              ``int`` or ``str``     See `Required Parameters`_
 ``nested_output``         ``bool``               See `Input/Output Parameters`_
========================  =====================  ========================================================================================


Learn/Apply Parameters
``````````````````````

General parameters related to Snekmer's learn and apply mode (``snekmer learn``, ``snekmer apply``) , wherein supervised models are trained via the workflow.

=============================  =====================  =========================================================================
     Parameter                    Type                 Description
=============================  =====================  =========================================================================
 ``save_apply_associations``     ``bool``              Save large optional output files containing all generated cosine similarity scores.
 ``conf_weight_modifier``        ``int``               Weighting modifer for updating confidence when adding data to an existing kmer count matrix.
 ``fragmentation``               ``bool``              Option to fragment training data with multiple sub-options listed below.
 ``version``                     ``str``               Choose 'absolute' or 'percent'. An absolute length of 50 would be 50 amino acids long.
 ``frag_length``                 ``int``               Length of fragment. Depending on "version", this is a percent or absolute length.
 ``min_length``                  ``int``               Minimum length of fragment that should be retained. Values less than this are discarded.
 ``location``                    ``str``               Choose 'start', 'end', or 'random'. This is where on a sequence a fragment is taken from.
 ``seed``                        ``int``               Choose any (random) seed for reproducible fragmentation.
 ``selection``                   ``str``               The method for selecting an annotation: 'top_hit', 'greatest_distance', or 'combined_distance'.
 ``threshold``                   ``str`` or ``None``             A family-specific threshold used for prediction filtering: None, 'Median', 'Mean', '90th Percentile', etc.
 ``weight_top``                  ``float``             When selection method is 'combined_distance', this is the weight for the top_hit method.
 ``weight_distance``             ``float``             When selection method is 'combined_distance', this is the weight for the greatest_distance method.
 ``apply_output``                ``str`` or ``None``   The output name for the apply results in single file format.
=============================  =====================  =========================================================================



