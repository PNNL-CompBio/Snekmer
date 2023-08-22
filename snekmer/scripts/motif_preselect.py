#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  2 10:32:01 2023

author: @tnitka
"""

import pickle
import snekmer as skm
import pandas as pd
import numpy as np
import gzip
import gc
from sklearn.svm import LinearSVC

# ---------------------------------------------------------
# Files and parameters
# ---------------------------------------------------------

config = snakemake.config

# ---------------------------------------------------------
# Run script
# ---------------------------------------------------------

# load input data
with open(snakemake.input.matrix, "rb") as f:
    data = pickle.load(f)

data.astype({'label': 'category'})
data.astype({'background': 'boolean'})

with gzip.open(snakemake.input.weights, "rb") as f:
    weights = pd.read_csv(f)
    
# prevent kmer NA being read as np.nan
if config["k"] == 2:
    weights["kmer"] = weights["kmer"].fillna("NA")


kmers = weights['kmer'].values    
scores = weights['sample'].values
family = skm.utils.get_family(
    skm.utils.split_file_ext(snakemake.input.weights)[0],
    regex=config["input_file_regex"],
)

del weights
gc.collect()
    
# set category label name (e.g. "family")
label = config["score"]["lname"] if str(config["score"]["lname"]) != "None" else "label"

# binary T/F for classification into family
family = skm.utils.get_family(snakemake.wildcards.nb, regex=config["input_file_regex"])
binary_labels = [True if value == family else False for value in data[label]]

# get and sort only unique labels
unique_labels = np.unique(data[label])
unique_labels.sort()
score_index = np.searchsorted(unique_labels, family)

svm = LinearSVC(class_weight="balanced", random_state=None, max_iter=1000000)
vecs=np.array(data["sequence_vector"].astype(str).str.strip('[]').str.split(",").tolist(), dtype='float')
svm.fit(vecs, data[label])

sequences = pd.DataFrame(vecs)

# del data, vecs
gc.collect()
scores = pd.DataFrame(svm.coef_)
unit_score = max(scores.iloc[score_index].values)
for i in range(len(scores.iloc[score_index].values)):
    scores.iloc[score_index, i] = scores.iloc[score_index, i]/unit_score

kmers = pd.Series(kmers)
while scores.iloc[score_index].lt(-0.2).sum()>0:
    # temp_scores = scores
    features = list()
    for i in range(len(scores.iloc[score_index].values)):
        if scores.iloc[score_index, i]<-0.1:
            features.append(i)
            
    scores.drop(scores.columns[features], axis=1, inplace=True)
    kmers.drop(features, inplace=True)
    kmers.index = np.arange(len(kmers.index))
    sequences.drop(sequences.columns[features], axis=1, inplace=True)
    sequences.columns = np.arange(len(sequences.columns))
    del features
    gc.collect()
            
    vecs = sequences.to_numpy()
    svm.fit(vecs, data[label])
    scores = pd.DataFrame(svm.coef_)
    unit_score = max(scores.iloc[score_index].values)
    for i in range(len(scores.iloc[score_index].values)):
        scores.iloc[score_index, i] = scores.iloc[score_index, i]/unit_score



del svm, data, vecs
gc.collect()
    
# save output
scores.iloc[score_index].to_csv(snakemake.output.data, index=False, compression="gzip")
kmers.to_csv(snakemake.output.kmers, index=False, compression="gzip")
sequences.to_csv(snakemake.output.vecs, index=False, compression="gzip")