# ---------------------------------------------------------
# Imports
# ---------------------------------------------------------
import pickle
from datetime import datetime

import numpy as np
import pandas as pd
import snekmer as skm
from sklearn.model_selection import train_test_split, StratifiedKFold

# ---------------------------------------------------------
# Files and Parameters
# ---------------------------------------------------------

config = snakemake.config

# ---------------------------------------------------------
# Run script
# ---------------------------------------------------------

# log script start time
start_time = datetime.now()
label = (
    config["score"]["lname"] if str(config["score"]["lname"]) != "None" else "label"
)  # e.g. "family"

# get kmers for this particular set of sequences
kmers = skm.io.load_pickle(snakemake.input.kmerobj)

# load vectorized seq data
kmerlist, df = skm.io.load_npz(snakemake.input.data)

# loading all other models and harmonize basis sets
# for i in range(len(data)):
#     df = data[i]
#     kmerlist = kmer_sets[i]
#     vecs = skm.utils.to_feature_matrix(df["sequence_vector"].values)
#     df["sequence_vector"] = kmers.harmonize(vecs, kmerlist).tolist()
#     data[i] = df


# create matrix of kmer counts
x, y = len(df["sequence_vector"]), len(df["sequence_vector"][0])
matrix = np.zeros(x * y).reshape((x, y))
for i in range(x):
    for j in range(y):
        value = df["sequence_vector"][i]
        value = value[j]
        matrix[i, j] = value

weight = config["score"]["background_weight"]
probas = skm.score.feature_class_probabilities(
    matrix.T, df["filename"], kmers=kmers.basis.basis
)
scores = weight * probas["score"].values

# save score weights
np.savez_compressed(
    output.scores,
    kmerlist=kmerlist,
    probas=probas,
    scores=scores,
)
