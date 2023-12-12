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

# create matrix of kmer counts
x, y = len(df["sequence_vector"]), len(df["sequence_vector"][0])
matrix = np.zeros(x * y).reshape((x, y))
for i in range(x):
    for j in range(y):
        value = df["sequence_vector"][i]
        value = value[j]
        matrix[i, j] = value

probas = skm.score.feature_class_probabilities(
    matrix.T, df["filename"], kmers=kmers.basis.basis
)

# save score weights
np.savez_compressed(
    snakemake.output.scores,
    kmerlist=kmerlist,
    probas=probas["probability"].values,
)
