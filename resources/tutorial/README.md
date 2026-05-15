# Snekmer Tutorial Notebooks

Run these notebooks from the `resources/tutorial/` directory.

| Notebook | Audience | Description |
|---|---|---|
| `snekmer_learn_apply_user_tutorial.ipynb` | **All users** | Runs `easy-learn-apply` via the command line and shows how to read, filter, and visualize results. Start here. |
| `snekmer_learn_apply_rules_reference.ipynb` | **Developers / advanced users** | Step-by-step walkthrough of every internal pipeline rule using Snekmer's Python APIs directly. For users who want to understand the method or integrate pipeline steps into their own code. |
| `snekmer_model_cluster_search_demo.ipynb` | **All users** | Demonstrates the Model, Cluster, and Search pipeline on a small demo dataset. |

## Setup

Activate your Snekmer virtual environment and register it as a Jupyter kernel (one-time):

```bash
source ~/snekmer_env/bin/activate
pip install ipykernel
python -m ipykernel install --user --name=snekmer
jupyter notebook
```

Select the **snekmer** kernel when opening any notebook.

## Rendered versions

Rendered (pre-run) versions of these notebooks are also available in the
[online documentation](https://snekmer.readthedocs.io/en/latest/tutorial/).
