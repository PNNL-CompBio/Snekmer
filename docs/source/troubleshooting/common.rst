Frequently Asked Questions
==========================

Some commonly encountered questions are addressed here. For more
detailed or specific questions that have not been included below, feel free to
`submit an issue on Github <https://github.com/PNNL-CompBio/Snekmer/issues>`_.

Installation Questions
----------------------

For errors encountered during the installation process,
unless installation is successful for all packages but
specifically fails during the installation of the Snekmer
package itself, we recommend consulting the `venv documentation <https://docs.python.org/3/library/venv.html>`_ or searching the
`Conda troubleshooting page <https://conda.io/projects/conda/en/latest/user-guide/troubleshooting.html>`_
for the issue, depending on which installation method you are using. If packages other than Snekmer are also failing
to install, or else the environment is not created successfully,
the installation issues likely involve either venv or the
user's individual configuration.

Intel Macs
``````````

For older Macs with processors using the Intel x86-64 architecture, it is **highly**
recommended to use Python 3.11–3.13 and run the following commands after creating
your virtual environment:

.. code-block:: bash

   python -m pip uninstall -y numba llvmlite
   python -m pip install --upgrade pip setuptools wheel
   python -m pip install --no-cache-dir --only-binary=:all: "llvmlite==0.44.0" "numba==0.61.0"
   python - <<'PY'
   import numba, llvmlite
   print("numba:", numba.__version__)
   print("llvmlite:", llvmlite.__version__)
   PY

.. _troubleshooting/ubuntu:

Ubuntu Users
````````````

Ubuntu users need to install additional system packages before installing Snekmer's
dependencies. Run these commands **after** creating your virtual environment:

.. code-block:: bash

   sudo apt-get install python3.12-venv
   sudo apt install gcc g++
   sudo apt install -y python3.12-dev build-essential
   pip install --upgrade pip setuptools wheel
   pip install Cython numpy
   pip install --no-build-isolation hdbscan

Troubleshooting Error Messages
------------------------------

If you encounter an error while using Snekmer, we recommend first
checking the `Snakemake FAQs page <https://snakemake.readthedocs.io/en/stable/project_info/faq.html>`_
for a solution. We have also listed some common error messages below.
If your error message cannot be solved, feel free to let us know via
`Github issues <https://github.com/PNNL-CompBio/Snekmer/issues>`_.

MissingInputException
`````````````````````

Non-Windows Systems
:::::::::::::::::::

Generally, this error type means that Snekmer is unable to detect input files,
and the input files may not follow the required directory structure.
For more information, including an example of the file structure
needed, see :ref:`getting_started-configuration`.

Windows Systems
:::::::::::::::

For Windows systems specifically, there is a `known issue <https://github.com/PNNL-CompBio/Snekmer/issues/60>`_
with handling unzipping files that will raise a ``MissingInputException``
and cause Snekmer to terminate under failure. We are aware of this issue
and are actively working on a resolution; in the meantime, we recommend
separately unzipping any zipped files prior to evaluation via Snekmer.

FileNotFoundError
`````````````````

This error usually indicates that the file system latency is not
sufficiently long. In order words, the system is taking longer to
create the file than Snakemake is allotting to wait for the file to
be created. To troubleshoot this issue, we recommend increasing
the ``--latency`` parameter (see :ref:`getting_started-all_options`).
The default is set to 30 seconds, but the parameter can be adjusted
to suit your individual system.

/bin/sh: line 0: cd: {PATH}: No such file or directory
``````````````````````````````````````````````````````

This is a `known Snakemake issue <https://github.com/snakemake/snakemake/issues/1546>`_ that
can occur when the wrong Snakemake version is installed. Snekmer requires
``snakemake==9.13``. Verify your version with ``snakemake --version`` and
reinstall the correct version if necessary:

.. code-block:: bash

   pip install "snakemake==9.13"

Error: Directory cannot be locked.
``````````````````````````````````

The full error message should provide further instructions, but this
error will appear when Snekmer has been unexpectedly terminated.
Run ``snekmer {mode} --unlock`` (note: this command will not execute the
workflow) before rerunning the workflow.

If the error persists, delete the ``.snakemake`` directory and try again.

AttributeError in _load_configfile
``````````````````````````````````
Typically, this error arises if the path to the config.yaml file is not
specified correctly. To resolve this error, check that your config.yaml
file is located in the same directory from which you are executing Snekmer.
You can also specify the location of the config.yaml file,
e.g. ``snekmer {mode} --configfile /path/to/config.yaml``, to fix the issue.

OSError: [Errno 86] Bad CPU type in executable
``````````````````````````````````````````````

This occurs when Snekmer tries to run using an incorrect scheduler. Resolve it by
passing ``--scheduler greedy``:

.. code-block:: bash

   snekmer learn --scheduler greedy --configfile=./config.yaml

General Usage Questions
-----------------------

My logs are very long, resulting in large log file sizes. How can I reduce this?
````````````````````````````````````````````````````````````````````````````````

By default, Snekmer logs all Snakemake output, including construction of the DAG
and information about individual jobs. However, the default settings will produce
very big log files if several files are being evaluated at once. To reduce the
verbosity of output logs, we recommend invoking the ``--quiet`` parameter
(see :ref:`getting_started-all_options`).

Snekmer model mode is not working.
``````````````````````````````````

If ``snekmer model`` is not building models as intended, check
the following:

1. At least two input files are in the input directory. Note
   that Snekmer will not run as intended without 2+ input files.
2. All configuration parameters have been correctly specified.
   For more details and parameter descriptions, refer to
   :ref:`config-main`. Verify that no misspellings or invalid
   parameter specifications have been entered.
3. The ``snekmer model`` command is executed from the directory
   containing the input directory, and that a config.yaml file
   has been placed in the top-level directory. Refer to
   :ref:`getting_started-configuration` for the file structure
   required for Snekmer.

Snekmer cluster mode is producing an unusual number of clusters.
````````````````````````````````````````````````````````````````

If Snekmer cluster results in an unexpected number of clusters,
we recommend tuning the parameter set used to generate the clusters.
Most likely, the parameters used to generate the clusters are too
generalized, or specific, for the given dataset. For instance, if
Snekmer determines only 1 cluster for a given protein sequence set of
many individual sequences, the parameters guiding the clustering
algorithm is likely not sensitive enough to differentiate the underlying
clusters. See :ref:`Parameter Selection <background-params>` for more details.


``easy-learn-apply`` Questions
------------------------------

The pipeline failed partway through. How do I re-run it?
`````````````````````````````````````````````````````````

If ``easy-learn-apply`` fails after the ``learn`` step but before ``apply`` completes,
Snakemake may leave a lock file or partial outputs in the apply workspace. The safest
recovery is to delete the output directory and re-run from scratch:

.. code-block:: bash

   rm -rf my_results/
   snekmer easy-learn-apply --train ... --query ... --ann ... --output my_results

Alternatively, if ``learn`` completed successfully you can re-run only the apply step
directly:

.. code-block:: bash

   snekmer apply -d my_results/apply

No sequences have a high-confidence prediction.
```````````````````````````````````````````````

Check the following:

1. **Score = 0 for most sequences**: this means the query sequences share no k-mers
   with any training family. Verify that the same ``--k`` and ``--alphabet`` were used
   for both train and query. If the query sequences are from a very different organism
   or are very divergent, consider lowering the ``--k`` value (e.g. ``--k 6``).

2. **Confidence is low but Score > 0**: the model may be undertrained. Try providing
   more training sequences per family (at least 20–50 is recommended).

3. **Verify your annotation file** matches your training FASTA headers. Run:

   .. code-block:: python

      import pandas as pd
      from Bio import SeqIO

      ann = pd.read_csv("annotations.ann", sep="\t")
      ids_in_ann = set(ann["id"])

      ids_in_fasta = {r.id.split("|")[1] for r in SeqIO.parse("training.fasta", "fasta")
                      if "|" in r.id}

      print("Matched:", len(ids_in_ann & ids_in_fasta))
      print("In ann but not FASTA:", len(ids_in_ann - ids_in_fasta))

I see an ``INPUT ERROR`` message before the pipeline starts.
``````````````````````````````````````````````````````````````

Snekmer validates your inputs before launching the pipeline. Common messages and fixes:

**Vocabulary too large**: e.g. ``alphabet '5' has 10 symbols, k=8 → 100,000,000 k-mers``:
The combination of alphabet and k-mer length would create more k-mers than can fit in memory.
Use a coarser alphabet (lower number) or a smaller ``--k`` value.

**frag-length less than k**: e.g. ``--frag-length (4) is less than k (8)``:
Fragments shorter than k produce no k-mers. Increase ``--frag-length`` to at least the value of ``--k``.

**Annotation file missing 'family' column**: the ``.ann`` file must be tab-separated with
exactly the columns ``id`` and ``family``. Check for comma separators or misspelled headers.

**No training sequences found in annotation file**: no sequence IDs in your FASTA files
matched any ``id`` in the ``.ann`` file. Verify that the accession extraction matches your
header format (Snekmer extracts the field between the first pair of ``|`` characters for
UniProt-style headers). If your headers lack ``|``, use ``--create-ann`` or build the
``.ann`` file manually using your actual header IDs.

I used ``--create-ann`` but got an error about no annotated sequences.
```````````````````````````````````````````````````````````````````````

``--create-ann`` requires training FASTA headers in the format
``>db|FAMILY_LABEL|seqid ...``; the family label must be the field
**between the first pair of** ``|`` **characters**. Headers that use a different
format (e.g. ``>seqid description``) cannot be parsed automatically.

To check your headers:

.. code-block:: bash

   head -1 training_sequences.fasta

If the header does not contain ``|`` characters, you will need to provide a
``.ann`` file using ``--ann`` instead. See the
:doc:`easy-learn-apply tutorial <../tutorial/snekmer_easy_learn_apply_tutorial>`
for the annotation file format.
