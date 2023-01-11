Installation
============

We recommend `Mamba <https://mamba.readthedocs.io/en/latest/installation.html>`_
for installation handling. `Conda <https://www.anaconda.com/download/>`_ can be
used as an alternative, but Conda can take a long time to resolve dependencies,
thus rendering installation via Conda
significantly slower than installation via Mamba. Mamba/Conda will
both manage dependencies and versioning, which simplifies the
process of installation.

If you already have Conda but wish to use Mamba for installation,
you can install Mamba by running the following:

.. code-block:: bash

  conda install -c conda-forge mamba

Install Snekmer via Mamba/Conda
-------------------------------

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

**Warning:** Installation of Snekmer using ``pip`` is not recommended due to the complexity
of dependencies associated with Snakemake. Mamba/Conda will handle these automatically,
whereas ``pip`` will not.

The ``pip`` implementation has not been fully tested, but users may attempt installation
using the included specifications:

.. code-block:: bash

  pip install -r requirements.txt
  pip install -e git+https://github.com/PNNL-CompBio/Snekmer#egg=snekmer

Install Snekmer via Docker
--------------------------

Snekmer has been installed into a  public docker image hosted on `Dockerhub <https://hub.docker.com/repository/docker/jjacobson95/snekmer_env>`_.
Usage requires the of installation of `Docker Desktop <https://docs.docker.com/desktop/>`_.
This container is intended to be used via an interactive shell. Here, we provide the simplest method of usage.

To download and run a container:

.. code-block:: bash

  docker pull jjacobson95/snekmer_env:v1.0
  docker run jjacobson95/snekmer_env:v1.0


To use the command line interface within the container:

.. code-block:: bash

  docker ps       # This will display <container ID>
  docker exec -it <container ID> /bin/bash


Additional ``docker`` commands can be used to copy data into the container or to mount it to a local directory.

The current version of the Docker image requires the ``--configfile`` flag
to be called explicitly (see :ref:`getting_started-all_options`).

**Note:** This container is designed to run indefinitely and should be stopped after use.

Notes on Using Snekmer Docker Image with Apple Silicon (M1/M2 Mac) Systems
``````````````````````````````````````````````````````````````````````````

The current Docker image of Snekmer has not been built for compatibility with ARM64
systems (e.g. Apple silicon systems with the Apple M1/M2 chip). To use the Docker image
on an Apple silicon system, Rosetta 2 is required. See the
`relevant Docker documentation page <https://docs.docker.com/desktop/install/mac-install/>`_
for more information.

Rosetta 2 can be installed on the command line using ``softwareupdate``:

.. code-block:: bash

  softwareupdate --install-rosetta


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
