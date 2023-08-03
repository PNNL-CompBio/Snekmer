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
import numpy as np
import gzip
import gc
# from typing import Any, Dict, List, Optional
# from sklearn.base import BaseEstimator, ClassifierMixin
# from sklearn.tree import DecisionTreeClassifier
# from sklearn.ensemble import RandomForestClassifier, AdaBoostClassifier
# from sklearn.linear_model import LogisticRegression  # LogisticRegressionCV
# from sklearn.model_selection import GridSearchCV, cross_validate
# from sklearn.pipeline import make_pipeline, Pipeline
from sklearn.svm import LinearSVC

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

# data=data[[label, 'sequence_vector', 'background']]

data.astype({'label': 'category'})
data.astype({'background': 'boolean'})

    
# with open(snakemake.input.kmerobj, "rb") as f:
#     kmers = pickle.load(f)
    
with gzip.open(snakemake.input.weights, "rb") as f:
    weights = pd.read_csv(f)
    
with gzip.open(snakemake.input.vecs, "rb") as f:
    vecs=pd.read_csv(f)
    vecs.to_numpy
    
# prevent kmer NA being read as np.nan
if config["k"] == 2:
    weights["kmer"] = weights["kmer"].fillna("NA")


with gzip.open(snakemake.input.kmers, "rb") as f:
    kmers = pd.read_csv(f)
    
scores = weights['sample'].values
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

# get and sort only unique labels
unique_labels = np.unique(data['label'])
unique_labels.sort()
score_index = np.searchsorted(unique_labels, family)

# get alphabet name
if config["alphabet"] in skm.alphabet.ALPHABET_ORDER.keys():
    alphabet_name = skm.alphabet.ALPHABET_ORDER[config["alphabet"]].capitalize()
else:
    alphabet_name = str(config["alphabet"]).capitalize()
  
# run permutations and score each
motif = skm.motif.SnekmerMotif()
svm = LinearSVC(class_weight="balanced", random_state=None, max_iter=1000000000)
# svm= skm.model.SnekmerModel(
#     model="svc",
#     model_params={
#         "kernel": "linear",
#         "class_weight": "balanced",
#         "decision_function_shape": "ovr",
#         "random_state": None,
#                   },
#     params={
#         "scaler": scorer.scaler,
#         "random_state": config["model"]["random_state"],
#         "solver": "liblinear",
#         "class_weight": "balanced",
#     },
# )
# perm_data = motif.permute(
motif.permute(    
    data, 
    family,
    # skm.utils.get_family(snakemake.wildcards.nb, regex=config["input_file_regex"]),
    label_col=label)
    
# del data
gc.collect()
# scorer.fit(
#     kmers,
#     # list(kmerobj.kmer_set.kmers),
#     # perm_data,
#     data,
#     family,
#     # skm.utils.get_family(snakemake.wildcards.nb, regex=config["input_file_regex"]),
#     label_col=label,
#     vec_col="sequence_vector",
#     **config["score"]["scaler_kwargs"],)
# vecs=np.array(perm_data["sequence_vector"].astype(str).str.strip('[]').str.split(",").tolist(), dtype='float')
# vecs=np.array(data["sequence_vector"].astype(str).str.strip('[]').str.split(",").tolist(), dtype='float')
# features_in=np.vstack((kmers, vecs))
# svm.fit(vecs, perm_data[label])
svm.fit(vecs, data[label])
    
del data, motif, vecs
gc.collect()
perm_scores = pd.DataFrame(svm.coef_)

del svm
gc.collect()

unit_score = max(perm_scores.iloc[score_index].values)
for i in range(len(perm_scores.iloc[score_index].values)):
    perm_scores.iloc[score_index, i] = perm_scores.iloc[score_index, i]/unit_score
    
# save output
perm_scores.iloc[score_index].to_csv(snakemake.output.data, index=False, compression="gzip")

# record script endtime
#skm.utils.log_runtime(snakemake.log[0], start_time)