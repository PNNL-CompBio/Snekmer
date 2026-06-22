Installation
============
We recommend installing Snekmer using Python's `venv <https://docs.python.org/3/library/venv.html>`_ module. Before installing Snekmer, check that you have Python 3.11 or later installed using the command:

.. code-block:: bash

	python --version

The next step is to create and activate a new virtual environment that will contain Snekmer and its dependencies. :ref:`Ubuntu users should install additional dependencies after creating the virtual environment. <troubleshooting/ubuntu>`

.. code-block:: bash

	python -m venv ~/snekmer_env
	source ~/snekmer_env/bin/activate

Once the virtual environment has been created and activated, clone and install Snekmer using the following commands.

.. code-block:: bash

	# Download from Github Main Branch
	git clone https://github.com/PNNL-CompBio/Snekmer.git
	# Alternatively, download a specific branch:
	# git clone -b pre_paper_updates https://github.com/PNNL-CompBio/Snekmer.git
	cd Snekmer
	pip install --no-cache -r requirements.txt
	pip install .

Verify the installation succeeded:

.. code-block:: bash

	snekmer --version

You should see the current version number printed (e.g. ``1.4.1``). If the command
is not found, ensure your virtual environment is activated.

Testing Installation
====================

A quick test of the installation can be performed by running the following commands.

To test Learn/Apply, run from the root of the Snekmer repository:

.. code-block:: bash

   snekmer easy \
       --train  resources/demo_sequences/learn_apply_inputs/learn \
       --query  resources/demo_sequences/learn_apply_inputs/apply/test_sequences_1.fasta \
       --ann    resources/demo_sequences/learn_apply_inputs/annotations/TIGRFAMs_annotation.ann \
       --output test_results

Alternatively, a scripted demo that runs ``snekmer learn`` then ``snekmer apply``
step-by-step is available:

.. code-block:: bash

   cd resources/learn_apply_demo
   python3 run_demo.py

To test Model/Cluster/Search:

.. code-block:: bash

   cd resources/model_cluster_search_demo
   python3 run_demo.py

Troubleshooting Notes
`````````````````````

The full version of Snakemake is
`incompatible with Windows <https://snakemake.readthedocs.io/en/stable/getting_started/installation.html#full-installation>`_.
On Windows, you will need to use ``mamba`` to create the environment from the
minimal Snakemake specification. If you do not have mamba installed, install it
first from `miniforge <https://github.com/conda-forge/miniforge>`_ (recommended)
or via ``conda install -n base -c conda-forge mamba``. Then:

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

Install Snekmer via Docker
--------------------------

Snekmer has been installed into a public docker image hosted on `Dockerhub <https://hub.docker.com/r/jjacobson95/snekmer>`_.
Usage requires the installation of `Docker Desktop <https://docs.docker.com/desktop/>`_.

To download the image from Dockerhub:

.. code-block:: bash

  docker pull jjacobson95/snekmer:latest


To use the command line interface within the container:

.. code-block:: bash

  docker run --rm -v "$(pwd)":/work jjacobson95/snekmer:latest {mode} {args}

The Docker image accepts the same modes (cluster, model, search, learn, apply, and easy) and command line arguments as the Snekmer command line interface.



Optional: Blazing Signature Filter (BSF)
-----------------------------------------

BSF is an optional performance dependency for ``snekmer cluster``. It is not
optional; Snekmer works without it. See :doc:`Advanced / Optional Dependencies <advanced>`
for installation instructions.
