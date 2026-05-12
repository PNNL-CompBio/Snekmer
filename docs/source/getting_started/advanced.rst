Advanced / Optional Dependencies
=================================

The dependencies described here are **not required** to run Snekmer.
They provide optional performance improvements for specific use cases.


.. _bsf-install:

Blazing Signature Filter (BSF) — faster clustering
----------------------------------------------------

The `Blazing Signature Filter <https://github.com/PNNL-CompBio/bsf-jaccard-py>`_
is an optional dependency used by ``snekmer cluster`` to compute pairwise Jaccard
distance matrices more efficiently. It applies to the following clustering methods
(set via the ``cluster.method`` config parameter):

- ``density-jaccard``
- ``hdensity-jaccard``
- ``agglomerative-jaccard`` *(default)*

If BSF is not installed, Snekmer automatically falls back to
``scipy.spatial.distance.pdist`` for Jaccard distance computation. Clustering
will produce identical results — BSF simply runs faster on large datasets.

**BSF is not compatible with Apple silicon (M1/M2/M3) systems.** See the
`known Apple silicon issues <https://github.com/PNNL-CompBio/Snekmer/issues/102>`_
for details.

Install GCC (required before installing BSF)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

BSF requires GCC 4.9 or later. Install it for your operating system:

**macOS**

.. code-block:: bash

   brew install gcc llvm libomp

After installing ``llvm``, Homebrew may print a "Caveats" message with additional
flags that need to be set. Follow those instructions to ensure GCC is correctly
resolved. A typical caveats message looks like:

.. code-block:: none

   If you need to have llvm first in your PATH, run:
     echo 'export PATH="/usr/local/opt/llvm/bin:$PATH"' >> ~/.zshrc

   For compilers to find llvm you may need to set:
     export LDFLAGS="-L/usr/local/opt/llvm/lib"
     export CPPFLAGS="-I/usr/local/opt/llvm/include"

**Windows / Linux / Unix**

See the
`BSF documentation <https://github.com/PNNL-CompBio/bsf-jaccard-py#install-gcc-49-or-newers>`_
for platform-specific GCC installation instructions.

Install BSF
~~~~~~~~~~~~

With GCC installed and the Snekmer virtual environment active:

.. code-block:: bash

   pip install git+https://github.com/PNNL-CompBio/bsf-jaccard-py#egg=bsf
