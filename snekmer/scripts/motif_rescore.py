"""
Created on Tue Apr 25 15:14:39 2023

author: @tnitka
"""

# ---------------------------------------------------------
# Imports
# ---------------------------------------------------------
import pickle
# from datetime import datetime

import snekmer as skm
import pandas as pd
# import numpy as np
import gzip
import gc
# from typing import Any, Dict, List, Optional
# from sklearn.base import BaseEstimator, ClassifierMixin
# from sklearn.tree import DecisionTreeClassifier
# from sklearn.ensemble import RandomForestClassifier, AdaBoostClassifier
# from sklearn.linear_model import LogisticRegression  # LogisticRegressionCV
# from sklearn.model_selection import GridSearchCV, cross_validate
# from sklearn.pipeline import make_pipeline, Pipeline
# from sklearn.svm import SVC

# ---------------------------------------------------------
# Files and parameters
# ---------------------------------------------------------

config = snakemake.config

# ---------------------------------------------------------
# Run script
# ---------------------------------------------------------

# log script start time
# start_time = datetime.now()

# with open(snakemake.log[0], "a") as f:
#     f.write(f"start time:\t{{start_time}}\n")
    
# load input data
with open(snakemake.input.matrix, "rb") as f:
    data = pickle.load(f)
    
# set category label name (e.g. "family")
label = config["score"]["lname"] if str(config["score"]["lname"]) != "None" else "label"

data.astype({'label': 'category'})
    
# with open(snakemake.input.kmers, "rb") as f:
#     kmers = f.readlines()
    
with gzip.open(snakemake.input.weights, "rb") as f:
    weights = pd.read_csv(f)
    
# prevent kmer NA being read as np.nan
if config["k"] == 2:
    weights["kmer"] = weights["kmer"].fillna("NA")

kmers = weights['kmer'].values    
# scores = weights['sample'].values
family = skm.utils.get_family(
    skm.utils.split_file_ext(snakemake.input.weights)[0],
    regex=config["input_file_regex"],
)
scorer = skm.score.KmerScorer()

del weights
gc.collect()

# set number of permutations to test
n_iter = (
    config["motif"]["n"]  
    )


# get kmers for this particular set of sequences
# with open(snakemake.input.kmerobj, "rb") as f:
#     kmerobj = pickle.load(f)

# binary T/F for classification into family
family = skm.utils.get_family(snakemake.wildcards.nb, regex=config["input_file_regex"])
binary_labels = [True if value == family else False for value in data[label]]

# get alphabet name
if config["alphabet"] in skm.alphabet.ALPHABET_ORDER.keys():
    alphabet_name = skm.alphabet.ALPHABET_ORDER[config["alphabet"]].capitalize()
else:
    alphabet_name = str(config["alphabet"]).capitalize()
  
# run permutations and score each
motif = skm.motif.SnekmerMotif()
perm_data = motif.permute(
    data, 
    family,
    # skm.utils.get_family(snakemake.wildcards.nb, regex=config["input_file_regex"]),
    label_col=label)
del data
gc.collect()
scorer.fit(
    kmers,
    # list(kmerobj.kmer_set.kmers),
    perm_data,
    family,
    # skm.utils.get_family(snakemake.wildcards.nb, regex=config["input_file_regex"]),
    label_col=label,
    vec_col="sequence_vector",
    **config["score"]["scaler_kwargs"],)
del perm_data
gc.collect()
perm_scores = pd.DataFrame((scorer.probabilities["sample"]))
# score_out = pd.DataFrame(perm_scores)

del scorer
gc.collect()
    
# save output
perm_scores.to_csv(snakemake.output.data, index=False, compression="gzip")

# record script endtime
#skm.utils.log_runtime(snakemake.log[0], start_time)