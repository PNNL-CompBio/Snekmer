Installation
============

We recommend using `Anaconda <https://www.anaconda.com/download/>`_
to handle installation. Anaconda will manage dependencies and
versioning, which simplifies the process of installation.

Install Snekmer
---------------

Use conda to install an environment from the YML file with Snekmer and
all of its dependencies. (Note: Users may either download the
`YML file <https://github.com/PNNL-CompBio/Snekmer/blob/main/environment.yml>`_
directly from the repository, or clone the repository beforehand
using the ``git clone`` command.)

.. code-block:: bash

	conda env create -f environment.yml


Note that if you want to use the optional Blazing Signature Filter (BSF) to
speed up clustering you must follow the BSF installation instructions below
and then you can use the alternate conda environment.

.. code-block:: bash

  conda env create -f environment_BSF.yml

After the install completes activate the conda environment

.. code-block:: bash

  conda activate snekmer

The package should now be ready to use!

Troubleshooting Notes
`````````````````````

If you are a Windows user and running into conflicts/errors when
creating the conda environment, you may need to install the minimal
version of Snakemake:

.. code-block:: bash

  conda create -n snekmer -c conda-forge -c bioconda -c numba python>=3.9 biopython matplotlib numpy>=1.22.3 numba>=0.56 scipy pandas seaborn snakemake-minimal==7.0 scikit-learn

.. Install Snekmer
.. ---------------

.. Activate the conda environment:

.. .. code-block:: bash

.. 	conda activate snekmer

.. Then, install Snekmer using pip (note: git clone step is optional
.. if you already have the repo cloned locally):

.. .. code-block:: bash

..   # option 1: clone repository (if you haven't already) and install
..   git clone https://github.com/PNNL-CompBio/Snekmer.git
..   pip install Snekmer

..   # option 2: direct install (no repository download required)
..   pip install git+https://github.com/PNNL-CompBio/Snekmer

(optional) Install GCC for BSF
------------------------------

The `Blazing Signature Filter <https://github.com/PNNL-Compbio/bsf-jaccard-py>`_
is a pairwise similarity algorithm that can optionally be used to efficiently
compute a distance matrix for Snekmer's clustering mode.

**Note that BSF is not required to run Snekmer.** For users that do not want
to use BSF for clustering, these instructions can be ignored.

In order for BSF to install correctly, GCC 4.9+ must be
installed on your system using the following instructions for the listed
operating systems. Once GCC is installed successfully, follow the remaining
setup steps.

Mac
```

Install GCC and the relevant dependencies using Homebrew.

.. code-block:: bash

  brew install gcc llvm libomp

After installing ``llvm``, some flags and your ``PATH`` variable may need to
be updated. Homebrew will output a "Caveats" message that may resemble the one
shown below:

.. code-block:: none

  To use the bundled libc++ please add the following LDFLAGS:
    LDFLAGS="-L/usr/local/opt/llvm/lib -Wl,-rpath,/usr/local/opt/llvm/lib"

  llvm is keg-only, which means it was not symlinked into /usr/local,
  because macOS already provides this software and installing another version in
  parallel can cause all kinds of trouble.

  If you need to have llvm first in your PATH, run:
    echo 'export PATH="/usr/local/opt/llvm/bin:$PATH"' >> ~/.zshrc

  For compilers to find llvm you may need to set:
    export LDFLAGS="-L/usr/local/opt/llvm/lib"
    export CPPFLAGS="-I/usr/local/opt/llvm/include"

You may follow these instructions to ensure GCC is correctly pulled as needed.

Windows or Linux/Unix
`````````````````````

Please refer to the
`BSF documentation <https://github.com/PNNL-Compbio/bsf-jaccard-py#install-gcc-49-or-newers>`_
for Linux/Unix or Windows instructions for installing GCC.

BSF Install for Snekmer Use
```````````````````````````
In the snekmer conda environment use the command

.. code-block:: bash

   pip install git+https://github.com/PNNL-Compbio/bsf-jaccard-py#egg=bsf
