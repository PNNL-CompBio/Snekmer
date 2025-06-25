.. _config-main:

Setting up User Configuration
=============================

To run Snekmer, the user must specify parameters in a configuration
file (.YAML). A template ``config.yaml`` file is included in the
`resources directory <https://github.com/PNNL-CompBio/Snekmer/tree/main/resources>`_.

The example YAML files included are:

* ``config.yaml``: Configuration file for running Snekmer
* ``clust.yaml``: (optional) Cluster configuration file for deploying Snekmer on a high-performance computing (HPC) cluster

Parameter Descriptions for ``config.yaml``
------------------------------------------

The base `config.yaml` file is required in order to run `snekmer model` or `snekmer cluster`.

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
 ``input_file_exts``       ``list``               File extensions to be considered as valid for input sequence files
 ``input_file_regex``      ``str`` or ``None``    Regular expression for parsing family/annotation identifiers from filenames
 ``nested_output``         ``bool``               If True, saves files into nested directory structure, i.e. `{save_dir}/{alphabet}/{k}`
========================  ====================  =========================================================================

Score Parameters
````````````````

General parameters related to how Snekmer calculates family scores for k-mers.

========================  =====================  =================================================================================
     Parameter                   Type             Description
========================  =====================  =================================================================================
 ``scaler_kwargs``         ``dict``               Keyword arguments to pass to k-mer scaler object
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
 ``method``                ``str``                Clustering method (options: ``"kmeans"``, ``"agglomerative"``,
                                                  ``"correlation"``, ``"density"``, ``"birch"``, ``"optics"``,
                                                  or ``"hdbscan"``)
 ``params``                ``dict``               Parameters to pass to the clustering algorithm
 ``cluster_plots``         ``bool``               If True, generates plots illustrating clustering results
 ``min_rep``               ``int`` or ``None``    Threshold for the minimum number of repetitions of a kmer within a set.
                                                  Kmers that do not meet this threshold are discarded.
 ``max_rep``               ``int`` or ``None``    Threshold for the maximum number of repetitions of a kmer within a set.
                                                  Kmers that do not meet this threshold are discarded.
 ``save_matrix``           ``bool``               If True, saves distance matrices (BSF). Not recommended for large datasets.
 ``dist_thresh``           ``int``                Distance threshold for BSF matrix
========================  ====================  ==============================================================================

Parameter Descriptions for ``clust.yaml``
-------------------------------------------

See `SLURM documentation <https://slurm.schedmd.com/sbatch.html>`_ for more information on cluster parameters.

Required Parameters for Snekmer Search
--------------------------------------

The following parameters are required in your config file for `snekmer search`.

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
````````````````

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
=============================  =====================  =========================================================================


Motif Parameters
````````````````
The following parameters are required for Snekmer's motif mode (``snekmer motif``), wherein feature selection is performed to find functionally relevant kmers.

========================  =====================  ==================================================================================
     Parameter                    Type            Description
========================  =====================  ==================================================================================
``n``                     ``int``                Number of label permutation and rescoring iterations to run for each input family.
========================  =====================  ==================================================================================
