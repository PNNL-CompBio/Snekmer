Quick Start
===========

Once Snekmer is installed, you can annotate sequences with a single command — no
configuration file or directory setup required.

Try it with the included demo data
------------------------------------

The repository includes demo data — 5,000 training proteins across 200 TIGRFAM families
and 3,000 query proteins — so you can run a complete example immediately after installing.
From the root of the Snekmer repository:

.. code-block:: bash

   snekmer easy-learn-apply \
       --train  resources/demo_sequences/learn_apply_inputs/learn \
       --query  resources/demo_sequences/learn_apply_inputs/apply/test_sequences_1.fasta \
       --ann    resources/demo_sequences/learn_apply_inputs/annotations/TIGRFAMs_annotation.ann \
       --output my_first_results

Results will be written to ``my_first_results/apply/snekmer_results.csv``.
See the :doc:`tutorial <../tutorial/snekmer_easy_learn_apply_tutorial>` for a
full walkthrough of this demo including output interpretation.

----

What you need
-------------

Three inputs:

1. **Training sequences**: a FASTA file or directory of FASTA files containing
   proteins with known family assignments.
2. **Query sequences**: a FASTA file or directory of FASTA files to annotate.
3. **Annotations**: one of:

   - ``--ann PATH``: a tab-separated file mapping sequence IDs to family labels:

     .. code-block:: text

        id              family
        A0A2D0MWR0      TIGR04183
        A0A1Y4R5C6      TIGR00722

     The ``id`` column must match the accession in your training FASTA headers.
     For UniProt-style headers (``>db|ACCESSION|name ...``), the accession is the
     field between the first pair of ``|`` characters.

     Family labels can be **any string** — database accessions (e.g. ``TIGR04183``),
     descriptive names (e.g. ``nitrogenase``), or numbers. Labels are case-sensitive.

   - ``--create-ann``: generate annotations automatically from training FASTA
     headers. Requires headers in the format ``>db|FAMILY_LABEL|seqid ...``, where
     the field between the first ``|`` pair is used as the family label.

Run
---

.. code-block:: bash

   snekmer easy-learn-apply \
       --train  path/to/training_sequences/ \
       --query  path/to/query_sequences.fasta \
       --ann    path/to/annotations.ann \
       --output results/

Or, if your FASTA headers encode family labels (``>db|FAMILY|seqid ...``):

.. code-block:: bash

   snekmer easy-learn-apply \
       --train  path/to/training_sequences/ \
       --query  path/to/query_sequences.fasta \
       --create-ann \
       --output results/

Snekmer will prompt for any missing inputs if flags are omitted.

Results
-------

The main output is a CSV file at ``results/apply/snekmer_results.csv``:

.. code-block:: text

   Sequence              Prediction   Score   delta   Confidence
   tr|A0A427BXE3|...     TIGR01023    0.198   0.08    1.000
   tr|A0A401ZGP4|...     TIGR00757    0.316   0.02    0.922
   ...

- **Prediction**: predicted family (highest cosine similarity to training profiles)
- **Score**: cosine similarity to the predicted family
- **delta**: gap between top and second-best score (larger = more certain)
- **Confidence**: calibrated probability the prediction is correct (0–1)

A confidence of **≥ 0.95** is a reliable starting threshold for high-quality annotations.

Next steps
----------

- :doc:`easy-learn-apply tutorial <../tutorial/snekmer_easy_learn_apply_tutorial>`:
  full walkthrough with demo data, output interpretation, and post-hoc evaluation.
- :doc:`Configuration <config>`: tune k-mer length, alphabet, scoring thresholds,
  and other parameters.
- :doc:`Usage <usage>`: full reference for all Snekmer modes.
