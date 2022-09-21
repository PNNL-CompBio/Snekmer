Frequently Asked Questions
==========================

Some commonly encountered questions are addressed here. For more
detailed or specific questions, feel free to `submit an issue on Github <https://github.com/PNNL-CompBio/Snekmer/issues>`_.

Installation
------------

Why is the ``conda install`` command is taking a long time?
```````````````````````````````````````````````````````````

Installation of ``snakemake`` via conda typically requires a fairly
long compile time. To avoid this, you can install Snakemake via
`Mamba <https://github.com/mamba-org/mamba>`_` (see the official
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


Usage
-----

Why I am encountering a ``MissingInputException`` when I try to run Snekmer?
````````````````````````````````````````````````````````````````````````````

Typically, this means that Snekmer is unable to detect input files,
and the input files may not follow the required directory structure.
For more information, including an example of the file structure
needed, see :ref:`getting_started-configuration`.

I cannot build models correctly.
````````````````````````````````

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