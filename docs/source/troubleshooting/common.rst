Frequently Asked Questions
==========================

Some commonly encountered questions are addressed here. For more
detailed or specific questions, feel free to `submit an issue on Github <https://github.com/PNNL-CompBio/Snekmer/issues>`_.

Installation
------------

The ``conda install`` command is taking a long time. Why?
`````````````````````````````````````````````````````````

Installation of ``snakemake`` via conda typically requires a fairly
long compile time. To avoid this, you can install Snakemake via
`Mamba <https://github.com/mamba-org/mamba>`_` (see the official
`Snakemake installation instructions <https://snakemake.readthedocs.io/en/stable/getting_started/installation.html>`_
for details).


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

*To be continued...*