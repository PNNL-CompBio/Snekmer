"""
Created on Tue Apr 25 15:14:39 2023

author: @tnitka
"""

# ---------------------------------------------------------
# Imports
# ---------------------------------------------------------
import pickle
from datetime import datetime

import snekmer as skm
import pandas as pd
import numpy as np
import gzip
from typing import Any, Dict, List, Optional
from sklearn.base import BaseEstimator, ClassifierMixin
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier, AdaBoostClassifier
from sklearn.linear_model import LogisticRegression  # LogisticRegressionCV
from sklearn.model_selection import GridSearchCV, cross_validate
from sklearn.pipeline import make_pipeline, Pipeline
from sklearn.svm import SVC

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
    
# with open(snakemake.input.kmers, "rb") as f:
#     kmers = f.readlines()
    
with gzip.open(snakemake.input.weights, "rb") as f:
    weights = pd.read_csv(f)

kmers = weights['kmer'].values    
scores = weights['sample'].values
family = skm.utils.get_family(
    skm.utils.split_file_ext(snakemake.input.weights)[0],
    regex=config["input_file_regex"],
)
scorer = skm.score.KmerScorer()

# set number of families
# n_families = (
#     config["motif"]["families"]
#     )

# set number of permutations to test
n_iter = (
    config["motif"]["n"]  
    )


# get kmers for this particular set of sequences
with open(snakemake.input.kmerobj, "rb") as f:
    kmerobj = pickle.load(f)
    
# set category label name (e.g. "family")
label = config["score"]["lname"] if str(config["score"]["lname"]) != "None" else "label"

# binary T/F for classification into family
family = skm.utils.get_family(snakemake.wildcards.nb)
binary_labels = [True if value == family else False for value in data[label]]

# prevent kmer NA being read as np.nan
if config["k"] == 2:
    scores["kmerobj"] = scores["kmerobj"].fillna("NA")
if config["k"] == 2:
    scores["kmers"] = scores["kmers"].fillna("NA")

# get alphabet name
if config["alphabet"] in skm.alphabet.ALPHABET_ORDER.keys():
    alphabet_name = skm.alphabet.ALPHABET_ORDER[config["alphabet"]].capitalize()
else:
    alphabet_name = str(config["alphabet"]).capitalize()

  
# run permutations and score each

score_matrix = pd.DataFrame({'kmer': kmers})
score_array = pd.DataFrame.to_numpy(score_matrix)
motif = skm.motif.SnekmerMotif()
for i in range(n_iter):
    perm_data = motif.permute(
        data, 
        skm.utils.get_family(snakemake.wildcards.nb, regex=config["input_file_regex"]),
        label_col=label)
    scorer.fit(
        list(kmerobj.kmer_set.kmers),
        perm_data,
        skm.utils.get_family(snakemake.wildcards.nb, regex=config["input_file_regex"]),
        label_col=label,
        vec_col="sequence_vector",
        **config["score"]["scaler_kwargs"],)
    perm_scores = np.reshape(pd.DataFrame.to_numpy(pd.DataFrame(scorer.probabilities["sample"])), (len(kmers), 1))

    score_array = np.hstack((score_array, perm_scores))
    
else:
    score_matrix=score_matrix.merge(
        pd.DataFrame(score_array), left_index=True, right_index=True
    )
    
output_matrix = motif.p_values(score_matrix, scores, n_iter)
    
# save output
# kmers.to_csv(snakemake.output.data, index=False, compression="gzip") # this should be redundant
output_matrix.to_csv(snakemake.output.p_values, index=False, compression="gzip")

# record script endtime
#skm.utils.log_runtime(snakemake.log[0], start_time)