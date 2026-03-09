========================================
Snekmer Learn/Apply Tutorial
========================================

This tutorial walks through the process of running the **Snekmer Learn** and
**Snekmer Apply** pipelines, from directory setup through vectorization and
annotation prediction. A detailed explanation of each step is available in ``resources/tutorial/snekmer_learn_apply_tutorial.ipynb``.

Overview
========

The basic process for running Snekmer Learn/Apply is as follows:

1. Verify that your file directory structure is correct and that the ``config.yaml`` file is in the top-level directory.

   - A ``config.yaml`` template is included in the Snekmer codebase at
     ``resources/learn_apply/config.yaml``.

2. Modify ``config.yaml`` with the desired parameters.

3. Use the command line to navigate to the directory containing both the
   ``config.yaml`` file and the ``input/`` directory.

4. Run ``snekmer learn``, then copy the appropriate outputs to a separate
   directory to run ``snekmer apply``.


Running the Snekmer Learn Pipeline
===================================

Configuration Parameters
~~~~~~~~~~~~~~~~~~~~~~~~~

Below is an example configuration dictionary.  These values typicaly live in ``config.yaml``, but can also be specified from the command line.

.. code-block:: python

   config = {
       # ...
       "save_apply_associations": False,

       # Fragmentation settings
       "fragmentation": False,
       "version": "absolute",
       "frag_length": 50,
       "min_length": 50,
       "location": "random",
       "seed": 999,

       # Confidence weighting
       "conf_weight_modifier": 20,

       # Selection methods: "top_hit", "greatest_distance", "combined_distance"
       "selection": "top_hit",

       # Threshold column from family_summary_stats.csv (or None)
       "threshold": "Median",

       # Weights (only relevant when using combined_distance)
       "weight_top": 0.5,
       "weight_distance": 0.5,

       # Output naming
       "apply_output": "snekmer_results.csv",
   }

Visualising Prediction Results
==============================

After running both **Learn** and **Apply**, the annotation predictions can be
visualised with a simple bar chart:

.. code-block:: python

   fig, ax = plt.subplots()
   # ... (build bar chart from results DataFrame) ...
   ax.set_title("Annotation Prediction Results")
   ax.set_ylabel("Sequences")
   h, l = ax.get_legend_handles_labels()
   u = dict(zip(l, h))
   ax.legend(u.values(), u.keys(), ncol=2, fontsize=8)
   plt.tight_layout()
   plt.show()
