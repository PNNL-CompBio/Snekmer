Setting up User Configuration
=============================

To run Snekmer, the user must specify parameters in a configuration
file (.YAML). A template ``config.yaml`` file is included in the
`resources directory <https://github.com/PNNL-CompBio/Snekmer/tree/main/resources>`_.

The example YAML files included are:

* ``config.yaml``: Configuration file for running Snekmer
* ``search.yaml``: Snekmer search mode configuration file; required for search mode, but not necessary if only model or cluster mode are desired
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
                                             (see `documentation <https://github.com/PNNL-CompBio/Snekmer/blob/main/snekmer/alphabet.py#L30>`_
                                             for more details). Alphabets may be specified by numbers 0-5 or by their names.
====================  ====================  ===================================================================================================

Alphabets 
:::::::::

Snekmer comes with 6 different recoding schemes, or "alphabets",
by which amino acid reduction (AAR) can occur. The available
options, including the option to skip recoding, are listed below.

=============  ===============  ======  ===========================================================================================
 Alphabet No.   Alphabet Name    Size                                         Description  
=============  ===============  ======  ===========================================================================================
 ``0``         ``hydro``          2      2-value hydrophobicity alphabet from Arnold, *et al.* 2009. PLoS Pathogens 5(4): e1000376
-------------  ---------------  ------  -------------------------------------------------------------------------------------------
 ``1``         ``standard``       7      "Standard" reduction alphabet from Arnold, *et al.* 2009. PLoS Pathogens 5(4): e1000376
-------------  ---------------  ------  -------------------------------------------------------------------------------------------
 ``2``         ``solvacc``        3      Solvent accessibility alphabet from Bacardit, et al. 2009. BMC Bioinformatics 10:6
-------------  ---------------  ------  -------------------------------------------------------------------------------------------
 ``3``         ``hydrocharge``    3      2-value hydrophobicity with charged residues as a third category; by @biodataganache
-------------  ---------------  ------  -------------------------------------------------------------------------------------------
 ``4``         ``hydrostruct``    3      2-value hydrophobicity with structural-breakers as a third category; by @biodataganache
-------------  ---------------  ------  -------------------------------------------------------------------------------------------
 ``5``         ``miqs``           10     MIQS alphabet; by @wichne
-------------  ---------------  ------  -------------------------------------------------------------------------------------------
 *n/a*         ``None``           20     No reduced alphabet
=============  ===============  ======  ===========================================================================================

Input/Output Parameters
```````````````````````

General parameters related to input and output sequences and/or files.

========================  ====================  =========================================================================
     Parameter                    Type            Description
========================  ====================  =========================================================================
 ``input_file_exts``       ``list``               File extensions to be considered as valid for input sequence files
 ``nested_output``         ``bool``               If True, saves files into nested directory structure, i.e. `{save_dir}/{alphabet}/{k}`
========================  ====================  =========================================================================

Score Parameters
````````````````

General parameters related to how Snekmer calculates family scores for k-mers.

========================  =====================  =================================================================================
     Parameter                   Type             Description
========================  =====================  =================================================================================
 ``scaler``                                       *currently inactive*
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
 ``n``                                            *currently inactive*
 ``cv``                    ``int``                Number of cross-validation folds for model evaluation
 ``use_score``                                    *currently inactive*
 ``random_state``          ``int`` or ``None``    Random state for model evaluation
========================  =====================  =========================================================================

Cluster Parameters
``````````````````

General parameters related to Snekmer's cluster mode (``snekmer cluster``), wherein unsupervised clusters are produced via the workflow.

========================  ====================  =========================================================================
     Parameter                    Type            Description
========================  ====================  =========================================================================
 ``method``                ``str``                Clustering method (options: ``"kmeans"``, ``"agglomerative"``,
                                                  ``"correlation"``, ``"density"``, ``"birch"``, or ``"optics"``)
 ``params``                ``dict``               Parameters to pass to the clustering algorithm
========================  ====================  =========================================================================
 

Parameter Descriptions for ``search.yaml``
------------------------------------------

The `search.yaml` file is required for `snekmer search`.

========================  =====================  ========================================================================================
     Parameter                     Type           Description
========================  =====================  ========================================================================================
 ``input_file_exts``       ``list``               File extensions to be considered as valid for input sequence files
 ``model_dir``             ``str``                Directory containing model object(s) (.model)
 ``basis_dir``             ``str``                Directory containing k-mer basis set(s) (.kmers)
 ``score_dir``             ``str``                Directory containing scoring object(s) (.scorer)
 ``k``                     ``int``                See `Required Parameters`_
 ``alphabet``              ``int`` or ``str``     See `Required Parameters`_
 ``nested_output``         ``bool``               See `Output Parameters`_
========================  =====================  ========================================================================================


Parameter Descriptions for ``clust.yaml``
-------------------------------------------

See `SLURM documentation <https://slurm.schedmd.com/sbatch.html>`_ for more information on cluster parameters.

