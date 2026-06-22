.. _Tutorial:
Tutorial
========

Learn/Apply Tutorials
---------------------

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - Tutorial
     - Description
   * - :doc:`Quick Start (easy) <snekmer_easy_learn_apply_tutorial>`
     - Recommended starting point. Single-command pipeline using the included demo data,
       with output interpretation and filtering guidance.
   * - :doc:`Full Pipeline Reference (learn / apply) <snekmer_learnapp_tutorial>`
     - For advanced use: manual workspace setup, adding new training data, reusing
       a trained model across multiple query sets.

Model / Cluster / Search Tutorial
----------------------------------

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - Tutorial
     - Description
   * - :doc:`Model / Cluster / Search <snekmer_demo>`
     - Supervised ML modeling, unsupervised clustering, and sequence search using the
       included demo data.


Evaluating Results
------------------

Snekmer writes an HTML summary report to ``output/Snekmer_<Mode>_Report.html``
for each completed run. Below is an example report from Apply:

.. image:: ../../../resources/images/apply_report_example.png
        :align: center
        :alt: Example Snekmer Apply html report

See the :ref:`Accessing Results <usage-results>` section for details on all output files.
