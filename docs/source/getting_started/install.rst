Installation
============
We recommend installing Snekmer using Python's `venv <https://docs.python.org/3/library/venv.html>`_ module. Before installing Snekmer, check that you have Python 3.9 or 3.10 installed using the command:

.. code-block:: bash

	python --version

The next step is to create and activate a new virtual environment that will contain Snekmer and its dependencies.

.. code-block:: bash

	python3 -m venv ~/snekmer_env
	source ~/snekmer_env/bin/activate

Once the virtual environment has been created and activated, clone and install Snekmer using the following commands.

.. code-block:: bash

	git clone https://github.com/PNNL-CompBio/Snekmer.git
	cd Snekmer
	pip install -r requirements.txt
	pip install .     # this installs local snekmer

Install Snekmer via Mamba/Conda
-------------------------------
It is also possible to use `Mamba <https://mamba.readthedocs.io/en/latest/installation.html>`_
for installation handling. `Conda <https://www.anaconda.com/download/>`_ can be
used as an alternative, but Conda can take a long time to resolve dependencies,
thus rendering installation via Conda
significantly slower than installation via Mamba. Mamba/Conda will
both manage dependencies and versioning, which simplifies the
process of installation.

The simplest method for installation is via the included YML file, which will create
a new environment containing Snekmer and all of its dependencies. Users may either
directly download the
`YML file <https://github.com/PNNL-CompBio/Snekmer/blob/main/environment.yml>`_
directly, or clone/fork the repository to obtain a local copy of the repository and all
included files.

.. code-block:: bash


	mamba env create -f environment.yml

Note that if you want to use the optional Blazing Signature Filter (BSF) to
speed up clustering you must follow the BSF installation instructions below
and then you can use the alternate conda environment.

.. code-block:: bash

  mamba env create -f environment_BSF.yml

After the install completes activate the conda environment

.. code-block:: bash

  conda activate snekmer

The package should now be ready to use!

Note that the instructions above can be replicated, substituting ``mamba``
for ``conda``, for users who wish to use Conda to manage installation.

Troubleshooting Notes
`````````````````````

The full version of Snakemake is
`incompatible with Windows <https://snakemake.readthedocs.io/en/stable/getting_started/installation.html#full-installation>`_.
Thus, you will need to install the environment specifications that
include only the minimal version of Snakemake:

.. code-block:: bash

  mamba env create -f environment_Windows.yml

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

Install Snekmer via pip
-----------------------

The ``pip`` implementation has not been fully tested, but users may attempt installation
directly from PyPI:

.. code-block:: bash
  pip install snekmer

Install Snekmer via Docker
--------------------------

Snekmer has been installed into a public docker image hosted on `Dockerhub <https://hub.docker.com/repository/docker/jjacobson95/snekmer_env>`_.
Usage requires the installation of `Docker Desktop <https://docs.docker.com/desktop/>`_.

To download the image from Dockerhub:

.. code-block:: bash

  docker pull jjacobson95/snekmer:latest


To use the command line interface within the container:

.. code-block:: bash

  docker run --rm   -v "$(pwd)":/data   -w /data   jjacobson95/snekmer:latest {mode} {args}    # Run Snekmer


The Docker image accepts the same modes (cluster, model, search, learn, apply, and motif) and command line arguments as the Snekmer command line interface.



(optional) Install GCC for BSF
------------------------------

The `Blazing Signature Filter <https://github.com/PNNL-CompBio/bsf-jaccard-py>`_
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

**Note:** BSF is not compatible with Apple silicon systems; see the ongoing log
of `known Apple silicon issues <https://github.com/PNNL-CompBio/Snekmer/issues/102>`_.

Windows or Linux/Unix
`````````````````````

Please refer to the
`BSF documentation <https://github.com/PNNL-CompBio/bsf-jaccard-py#install-gcc-49-or-newers>`_
for Linux/Unix or Windows instructions for installing GCC.

BSF Install for Snekmer Use
```````````````````````````
In the snekmer conda environment use the command

.. code-block:: bash

   pip install git+https://github.com/PNNL-CompBio/bsf-jaccard-py#egg=bsf
