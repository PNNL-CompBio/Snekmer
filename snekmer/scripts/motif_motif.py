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
start_time = datetime.now()
label = (
    config["model"]["lname"] if str(config["model"]["lname"]) != "None" else "Label"
    )
with open(snakemake.log[0], "a") as f:
    f.write(f"start time:\t{start_time}\n")
    
# load input data
with open(snakemake.input.matrix, "rb") as f:
    data = pickle.load(f)
    
kmers = skm.io.load_pickle(snakemake.input.kmerobj)
weights = pd.read_csv(snakemake.input.weights)
scores = weights[0, :] #TODO check whether this is the correct column
family = skm.utils.get_family(
    skm.utils.split_file_ext(snakemake.input.weights)[0],
    regex=config["input_file_regex"],
)
scorer = skm.score.KmerScorer()

#set number of permutations to test
n_iter = (
    config["motif"]["n"]  
    )


# get kmers for this particular set of sequences
with open(snakemake.input.kmerobj, "rb") as f:
    kmer = pickle.load(f) 
    
# set category label name (e.g. "family")
label = config["score"]["lname"] if str(config["score"]["lname"]) != "None" else "label"

# prevent kmer NA being read as np.nan
if config["k"] == 2:
    scores["kmer"] = scores["kmer"].fillna("NA")

# get alphabet name
if config["alphabet"] in skm.alphabet.ALPHABET_ORDER.keys():
    alphabet_name = skm.alphabet.ALPHABET_ORDER[config["alphabet"]].capitalize()
else:
    alphabet_name = str(config["alphabet"]).capitalize()

  
# run permutations and score each
input_matrix = snakemake.input.matrix
score_matrix = np.reshape(np.array(kmers), (len(kmers),1))
labels = input_matrix[:, 0] # TODO check whether this is the correct column
for i in range(n_iter):
    perm_data = skm.motif.SnekmerMotif.permute(input_matrix, labels)
    scorer.fit(kmers, perm_data, labels)
    perm_scores = scorer.scores
    score_matrix = np.append(score_matrix, perm_scores, 1)
    
# save output
kmers.to_csv(snakemake.output.data, index=False, compression="gzip")
score_matrix.to_csv(snakemake.output.p_values, index=Fals, compression="gzip")

# record script endtime
skm.utils.log_runtime(snakemake.log[0], start_time)