Installation
============

We recommend using `Anaconda <https://www.anaconda.com/download/>`_
to handle installation. Anaconda will manage dependencies and
versioning, which simplifies the process of installation.

Setup
-----

Use conda to install an environment from the YML file with all
required dependencies. (Note: Users may either download the
`YML file <https://github.com/PNNL-CompBio/Snekmer/blob/main/environment.yml>`_
directly from the repository, or clone the repository beforehand
using the ``git clone`` command.)

.. code-block:: bash

	conda env create -f environment.yml

Troubleshooting Notes
`````````````````````

If you are a Windows user and running into conflicts/errors when
creating the conda environment, you may need to install the minimal
version of Snakemake:

.. code-block:: bash
  
  conda create -n snekmer -c conda-forge -c bioconda biopython matplotlib numpy pandas seaborn snakemake-minimal scikit-learn

Install Snekmer
---------------

Activate the conda environment:

.. code-block:: bash

	conda activate snekmer

Then, install Snekmer using pip (note: git clone step is optional
if you already have the repo cloned locally):

.. code-block:: bash

  # option 1: clone repository (if you haven't already) and install
  git clone https://github.com/PNNL-CompBio/Snekmer.git
  pip install Snekmer

  # option 2: direct install (no repository download required)
  pip install git+https://github.com/PNNL-CompBio/Snekmer

The package should now be ready to use!

