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
from sklearn.svm import LinearSVC
# import glob
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
    
# with open(snakemake.input.kmers, "rb") as f:
#     kmers = f.readlines()
    
with gzip.open(snakemake.input.weights, "rb") as f:
    weights = pd.read_csv(f)
    
# set category label name (e.g. "family")
label = config["score"]["lname"] if str(config["score"]["lname"]) != "None" else "label"
    
# with open(snakemake.input.model, "rb") as f:
#     model = pickle.load(f)

svm = LinearSVC(class_weight="balanced", random_state=None, max_iter=1000000)
vecs=np.array(data["sequence_vector"].astype(str).str.strip('[]').str.split(",").tolist(), dtype='float')
svm.fit(vecs, data[label])
    
# prevent kmer NA being read as np.nan
if config["k"] == 2:
    weights["kmer"] = weights["kmer"].fillna("NA")

kmers = weights['kmer'].values    
coeffs = pd.DataFrame(svm.coef_)
# scores = weights['sample'].values
family = skm.utils.get_family(
    skm.utils.split_file_ext(snakemake.input.weights)[0],
    regex=config["input_file_regex"],
)
unique_labels = np.unique(data['label'])
unique_labels.sort()
score_index = np.searchsorted(unique_labels, family)
scores = coeffs.iloc[score_index]
scorer = skm.score.KmerScorer()

del svm
gc.collect()

unit_score = max(scores)
for i in range(len(scores)):
    scores.iloc[i] = scores.iloc[i]/unit_score

# set number of permutations to test
n_iter = (
    config["motif"]["n"]  
    )


# get kmers for this particular set of sequences
# with open(snakemake.input.kmerobj, "rb") as f:
#     kmerobj = pickle.load(f)


# binary T/F for classification into family
family = skm.utils.get_family(snakemake.wildcards.nb)
binary_labels = [True if value == family else False for value in data[label]]

# get alphabet name
if config["alphabet"] in skm.alphabet.ALPHABET_ORDER.keys():
    alphabet_name = skm.alphabet.ALPHABET_ORDER[config["alphabet"]].capitalize()
else:
    alphabet_name = str(config["alphabet"]).capitalize()

  
# run permutations and score each

score_matrix = pd.DataFrame({'kmer': kmers})
score_array = pd.DataFrame.to_numpy(score_matrix)
motif = skm.motif.SnekmerMotif()
for file in snakemake.input.perm_scores:
    with gzip.open(file) as f:
        perm_scores = pd.DataFrame.to_numpy(pd.read_csv(f))
    score_array = np.hstack((score_array, perm_scores))
    
else:
    score_array = np.delete(score_array, 0, 1)
    score_matrix=score_matrix.merge(
        pd.DataFrame(score_array), left_index=True, right_index=True
    )
    
output_matrix = motif.p_values(score_matrix, scores, n_iter)
output_matrix = output_matrix.astype({'kmer': 'str', 'real score': 'float32', 'false positives': 'int32', 'n': 'int32', 'p': 'float32'})
output_matrix.sort_values(by=['p', 'real score'], ascending=[True, False], inplace=True)
    
# save output
score_matrix.to_csv(snakemake.output.data, index=False, compression="gzip")
output_matrix.to_csv(snakemake.output.p_values, index=False, compression="gzip")

# record script endtime
#skm.utils.log_runtime(snakemake.log[0], start_time)