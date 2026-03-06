# Snekmer: A scalable pipeline for protein sequence fingerprinting using amino acid recoding (AAR)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7662597.svg)](https://doi.org/10.5281/zenodo.7662597)
[![CI](https://github.com/PNNL-CompBio/Snekmer/actions/workflows/action.yml/badge.svg)](https://github.com/PNNL-CompBio/Snekmer/actions)
[![Documentation Status](https://readthedocs.org/projects/snekmer/badge/?version=latest)](https://snekmer.readthedocs.io/en/latest/?badge=latest)
[![Snakemake](https://img.shields.io/badge/snakemake-=7.0.0-brightgreen.svg?style=flat)](https://snakemake.readthedocs.io)

Snekmer is a software package designed to reduce the representation of protein sequences
by combining amino acid reduction (AAR) with the kmer approach. Based on the AAR-kmer representations,
Snekmer subsequently (1) clusters sequences using various unsupervised clustering algorithms,
(2) generates supervised machine learning models, or (3) searches sequences against pre-trained models
to determine probabilistic annotations.

<p align="center">
  <img align="center" src="resources/images/snekmer_workflow.svg">
</p>

There are six operation modes for Snekmer: `cluster`, `model`, `search`, `learn`, `apply`, and `motif`.

---
**Clustering, Modeling, and Search Workflows** 

**Cluster mode:** The user supplies files containing sequences in an appropriate format (e.g. FASTA).
Snekmer applies the relevant workflow steps and outputs the resulting clustering results in tabular form (.CSV),
as well as the cluster object itself (.cluster). Figures are also generated (e.g. t-SNE, UMAP) to help the user
contextualize their results.

**Model mode:** The user supplies files containing sequences in an appropriate format (e.g. FASTA).
Snekmer applies the relevant workflow steps and outputs the resulting models as objects (.model). Snekmer also
displays K-fold cross validation results in the form of figures (AUC ROC and PR AUC curves) and a table (.CSV).

**Search mode:** The user supplies files containing sequences in an appropriate format (e.g. FASTA)
and the models they wish to search their sequences against. Snekmer applies the relevant workflow steps
and outputs a table for each file containing model annotation probabilities for the given sequences.
   
---
**Learn/Apply Annotation Workflows**

**Learn mode:** The user supplies files containing sequences in an appropriate format (e.g. FASTA)
as well as an annotation file. Snekmer generates a kmer counts matrix with the summed kmer distribution
of each annotation recognized from the sequence ID. Snekmer then performs a self-evaluation to assess
confidence levels. There are two outputs, a counts matrix, and a global confidence distribution.

**Apply mode:** The user supplies files containing sequences in an appropriate format (e.g. FASTA)
and the outputs received from Learn. Snekmer uses cosine distance to predict the annotation of each
sequence from the kmer counts matrix. The output is a table for each file containing sequence annotation
predictions with confidence levels.

---
**In Development**

**Motif mode:** The user supplies files containing sequences in an appropriate format (e.g. FASTA)
and the outputs received from Model. Snekmer performs a feature selection workflow to produce a
list of motifs ordered by degree of conservation and a classification model using the selected features (.model).

--------------------------------------------------------------------------------
Installation
--------------------------------------------------------------------------------

Requirements:
- Python 3.11+
- A working C/C++ toolchain may be required for some dependencies depending on your platform.
- Snekmer orchestrates workflows via Snakemake (installed as a Python dependency).

Option A: Install from source (recommended for development)
1) Create and activate a virtual environment:
```
  python -m venv ~/snekmer_env
  source ~/snekmer_env/bin/activate     # bash/zsh
  # or:
  source ~/snekmer_env/bin/activate.csh # csh
```
2) Clone and install:
```
  git clone https://github.com/PNNL-CompBio/Snekmer.git
  cd Snekmer
  pip install -r requirements.txt
  pip install .
```
Verify installation:
```
  snekmer -h
```
--------------------------------------------------------------------------------
Quick Start: Running Snekmer
--------------------------------------------------------------------------------

Snekmer workflows are run via the `snekmer` CLI. Most users will run using a YAML configuration file
(`config.yaml`). Please view the documentation site for full instructions as these brief examples are not
sufficient to run all workflows in all environments.
Working directory and required inputs (concise)  

Most runs use:
```
  input/            sequence files (FASTA/FAA/FNA/etc.)
  config.yaml       workflow parameters (start from resources/config.yaml)
```
Learn/Apply also uses:
```
  annotations/{file}.ann
```

Requiremented Inputs
----------------------------------------

**Cluster**  
  ```
    - input/
    - config.yaml
```
  
**Model**  
```
    - input/
    - config.yaml
```
  
**Search**  
```
    - input/
    - config.yaml
    - model_dir/*.model
    - basis_dir/*.kmers
    - score_dir/*.scorer
```
  
**Learn**  
```
    - input/                         (training FASTAs)
    - annotations/{file}.ann
    - config.yaml
```
  
**Apply**  
```
    - input/                         (query FASTAs)
    - config.yaml
    - counts/kmer_counts_total.csv                 # From Learn
    - confidence/global_confidence_scores.csv      # From Learn
    - stats/family_summary_stats.csv               # From Learn
```

Basic run commands
------------------
```
  snekmer {mode} --cores 2 --configfile ./config.yaml
```
Help:
```
  snekmer -h
  snekmer {mode} -h
```

Docs:
  https://snekmer.readthedocs.io
   
--------------------------------------------------------------------------------
Demos and Tutorials
--------------------------------------------------------------------------------

Snekmer includes minimal demos and Jupyter tutorials under the `resources/` directory.

Two end-to-end demos (each has a single entrypoint script named `run_demo.py`):
1) Model -> Cluster -> Search demo:
   `resources/model_cluster_search_demo/run_demo.py`

2) Learn -> Apply demo:
   `resources/learn_apply_demo/run_demo.py`

Run a demo (example):
```
  cd resources/model_cluster_search_demo
  python run_demo.py
```

Tutorial notebooks:
- `resources/tutorial/snekmer_demo.ipynb`
- `resources/tutorial/snekmer_learn_apply_tutorial.ipynb`

--------------------------------------------------------------------------------
Citation Guidance
--------------------------------------------------------------------------------

1. McDermott, Jason E., Chang, Christine H., Jerger, Abby, Nelson, William B., & Jacobson, Jeremy R. (2023).
   Snekmer: A scalable pipeline for protein sequence fingerprinting using amino acid recoding (AAR) (v1.0.3).
   Zenodo. https://doi.org/10.5281/zenodo.7662597

2. Christine H Chang, William C Nelson, Abby Jerger, Aaron T Wright, Robert G Egbert, Jason E McDermott,
   Snekmer: a scalable pipeline for protein sequence fingerprinting based on amino acid recoding,
   Bioinformatics Advances, Volume 3, Issue 1, 2023, vbad005. https://doi.org/10.1093/bioadv/vbad005

--------------------------------------------------------------------------------
Maintainers
--------------------------------------------------------------------------------

Snekmer was written and is maintained by the following PNNL development team:
Christine Chang, Jeremy Jacobson, Abby Jerger, Tara Nitka, Bill Nelson, and Jason McDermott.

--------------------------------------------------------------------------------
License
--------------------------------------------------------------------------------

BSD 3-Clause License

Copyright (c) 2021, Pacific Northwest National Laboratory All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


--------------------------------------------------------------------------------
Disclaimers
--------------------------------------------------------------------------------

```
This material was prepared as an account of work sponsored by an agency of the United States Government. Neither the United States Government nor the United States Department of Energy, nor Battelle, nor any of their employees, nor any jurisdiction or organization that has cooperated in the development of these materials, makes any warranty, express or implied, or assumes any legal liability or responsibility for the accuracy, completeness, or usefulness or any information, apparatus, product, software, or process disclosed, or represents that its use would not infringe privately owned rights.

Reference herein to any specific commercial product, process, or service by trade name, trademark, manufacturer, or otherwise does not necessarily constitute or imply its endorsement, recommendation, or favoring by the United States Government or any agency thereof, or Battelle Memorial Institute. The views and opinions of authors expressed herein do not necessarily state or reflect those of the United States Government or any agency thereof.

PACIFIC NORTHWEST NATIONAL LABORATORY operated by BATTELLE for the UNITED STATES DEPARTMENT OF ENERGY under Contract DE-AC05-76RL01830
```
