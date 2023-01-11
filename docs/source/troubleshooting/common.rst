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
package itself, we recommend searching the
`Conda troubleshooting page <https://conda.io/projects/conda/en/latest/user-guide/troubleshooting.html>`_
for the issue. If packages other than Snekmer are also failing
to install, or else the environment is not created successfully,
the installation issues likely involve either Conda or the
user's individual configuration.


Why is the ``conda install`` command is taking a long time?
```````````````````````````````````````````````````````````

Installation of ``snakemake`` via conda typically requires a fairly
long compile time. To avoid this, you can install Snakemake via
`Mamba <https://github.com/mamba-org/mamba>`_ (see the official
`Snakemake installation instructions <https://snakemake.readthedocs.io/en/stable/getting_started/installation.html>`_
for details).

Why am I getting a ``CondaEnvException: Pip failed``? Does this mean my Snekmer installation failed?
````````````````````````````````````````````````````````````````````````````````````````````````````

This error appears when `BSF <https://github.com/PNNL-CompBio/bsf-jaccard-py>`_
is unable to be installed on a user's system. However, conda still
installs Snekmer and its dependencies without BSF, despite the error.
Since BSF is optional for Snekmer, Snekmer should still run without
any issues.

To verify this, run ``conda activate <ENV_NAME>`` (i.e. ``conda activate snekmer``),
followed by ``conda list``, to check that the environment was
created successfully and that Snekmer was indeed installed.


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

This is a `known issue with Snakemake v7.3.0+ <https://github.com/snakemake/snakemake/issues/1546>`_.
Check your Snakemake version and reinstall a lower version if necessary
(we recommend Snakemake v7.0).

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