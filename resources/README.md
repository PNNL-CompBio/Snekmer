# Resources

The `resources/` directory contains example configuration files, runnable demos, and Jupyter notebook tutorials to help you validate an installation.  
  
--------------------------------------------------------------------------------
Configuration files (YAML)  
--------------------------------------------------------------------------------

Example YAML files:  
   
- `resources/config.yaml`    
  Default configuration template for running Snekmer workflows.  
  
- `resources/clust.yaml` (optional)  
  Example cluster configuration for running Snakemake on an HPC scheduler (e.g., SLURM).   

  Notes:
  - `partition` and `account` should not be left blank in `clust.yaml`.  
  - `ntasks` is used as `--ntasks-per-node` in the SLURM submit command.  
  - This file is only needed for SLURM/HPC runs (see "Running on SLURM" below).  

These YAML files are useful if you prefer a config enabled (reccomended) workflow or want to run on HPC.  

--------------------------------------------------------------------------------
Running on SLURM (optional)  
--------------------------------------------------------------------------------

Snekmer can submit jobs to SLURM via Snakemake when you provide a cluster config.  
  
General pattern:  
```
  snekmer <mode> --clust resources/clust.yaml --jobs <N> [other arguments]
```

Example (cluster mode):  
```
  snekmer cluster --clust resources/clust.yaml --jobs 10
```

Key flags:  
```
--clust PATH  
  Path to a cluster configuration YAML (e.g., resources/clust.yaml).  
  
--jobs N  
  Maximum number of simultaneous jobs Snakemake will submit.  
  
--cores N  
  Maximum cores for scheduling/parallelism (behavior depends on Snakemake execution mode and rules).  
```

--------------------------------------------------------------------------------
Demos (end-to-end examples)  
--------------------------------------------------------------------------------

Snekmer includes two end-to-end demos. Each demo is a self-contained workspace with a
single entrypoint script named `run_demo.py`.  
  
Prerequisite:
- Install Snekmer and ensure `snekmer` is on your PATH (e.g., `pip install .` from the repo root).  
  
Demo 1: Model -> Cluster -> Search
----------------------------------

**What it does:**  
- Resets the demo workspace  
- Copies example protein FASTA inputs into the demo `input/` directory  
- Runs:  
  - `snekmer cluster`  
  - `snekmer model`  
  - `snekmer search`  
- Optionally gathers model artifacts into `output/example-model/` for convenience  
  
**Location:**  
- Demo workspace: `resources/model_cluster_search_demo/`  
- Inputs: `resources/demo_sequences/model_cluster_search_inputs/`  

**Run Locally:**  
```  
  cd resources/model_cluster_search_demo  
  python run_demo.py  
```  

**Outputs:**  
Artifacts are written under:  
- `resources/model_cluster_search_demo/output/`    
  Common subdirectories include: `cluster/`, `model/`, `scoring/`, `search/`, `vector/`, `kmerize/`  
  HTML reports may also be generated depending on your workflow settings:  
    - `Snekmer_Cluster_Report.html`  
    - `Snekmer_Model_Report.html`  
    - `Snekmer_Search_Report.html`  
  Example Output folder:  
    - `output/example-model/`  
  
Demo 2: Learn -> Apply
----------------------
  
**What it does:**  
- Sets up two workspaces:  
  - `learn/`  (training)  
  - `apply/`  (prediction)  
- Copies training FASTAs + annotation file into the learn workspace  
- Runs:  
  - `snekmer learn`  
- Copies required handoff files from learn -> apply  
- Copies a test FASTA into the apply workspace  
- Runs:  
  - `snekmer apply`  
  
**Location:**  
- Demo workspace: `resources/learn_apply_demo/`  
- Inputs: `resources/demo_sequences/learn_apply_inputs/`  
  - `annotations/`  
  - `learn/`   (training FASTAs)  
  - `apply/`   (test FASTA)  
  
**Run Locally:**  
```
  cd resources/learn_apply_demo  
  python run_demo.py  
```
  
**Outputs:**  
Artifacts are written under:  
- `resources/learn_apply_demo/learn/output/`  
- `resources/learn_apply_demo/apply/output/`  
A single-file apply results table may be created (depending on your workflow settings):  
- `resources/learn_apply_demo/apply/snekmer_results.csv`  
  
--------------------------------------------------------------------------------
Tutorials (Jupyter notebooks)
--------------------------------------------------------------------------------

The `resources/tutorial/` folder contains Jupyter notebooks that walk through Snekmer
concepts and workflows with narrative explanations.  
  
**Tutorial 1: Snekmer overview**  
- Notebook: `resources/tutorial/snekmer_demo.ipynb`  
- Purpose: introduces the overall workflow and how Snekmer model/cluster/search modes fit together  
  
**Run Locally:**    
```
  cd resources/tutorial  
  jupyter notebook
  # open snekmer_demo.ipynb
```
  
**Tutorial 2: Learn -> Apply workflow**  
- Notebook: resources/tutorial/snekmer_learn_apply_tutorial.ipynb  
- Purpose: explains Learn/Apply mode, expected inputs/outputs, and interpretation of results  
  
**Run Locally:** 
```
  cd resources/tutorial  
  jupyter notebook  
  # open snekmer_learn_apply_tutorial.ipynb  
```

--------------------------------------------------------------------------------
Troubleshooting demo runs  
--------------------------------------------------------------------------------
  
If a demo fails:
  
1) Confirm the CLI is installed and visible:
```
  snekmer -h  
  snekmer cluster -h
```
2) Re-run with additional debug output:
```
  snekmer {mode} --verbose  
```
3) Clear the demo workspace and retry:  
  Each run_demo.py resets its demo directories, but you can also remove `output/` and  
  `.snakemake/` manually if needed.  
  
For more troubleshooting guidance, see the full [documentation](https://snekmer.readthedocs.io) site.
