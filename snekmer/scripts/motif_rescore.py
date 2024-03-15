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

# set category label name (e.g. "family")
label = config["score"]["lname"] if str(config["score"]["lname"]) != "None" else "label"

data.astype({"label": "category"})
data.astype({"background": "boolean"})

with gzip.open(snakemake.input.weights, "rb") as f:
    weights = pd.read_csv(f)

with gzip.open(snakemake.input.vecs, "rb") as f:
    vecs = pd.read_csv(f)
    vecs.to_numpy

# prevent kmer NA being read as np.nan
if config["k"] == 2:
    weights["kmer"] = weights["kmer"].fillna("NA")


with gzip.open(snakemake.input.kmers, "rb") as f:
    kmers = pd.read_csv(f)

scores = weights["sample"].values
family = skm.utils.get_family(
    skm.utils.split_file_ext(snakemake.input.weights)[0],
    regex=config["input_file_regex"],
)
scorer = skm.score.KmerScorer()

del weights
gc.collect()

# set number of permutations to test
n_iter = config["motif"]["n"]

# binary T/F for classification into family
family = skm.utils.get_family(snakemake.wildcards.nb, regex=config["input_file_regex"])
binary_labels = [True if value == family else False for value in data[label]]

# get and sort only unique labels
unique_labels = np.unique(data["label"])
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
motif.permute(data, family, label_col=label)

# del data
gc.collect()
svm.fit(vecs, data[label])

del data, motif, vecs
gc.collect()
perm_scores = pd.DataFrame(svm.coef_)

del svm
gc.collect()

# normalize kmer weights
unit_score = max(perm_scores.iloc[score_index].values)
for i in range(len(perm_scores.iloc[score_index].values)):
    perm_scores.iloc[score_index, i] = perm_scores.iloc[score_index, i] / unit_score

# save output
perm_scores.iloc[score_index].to_csv(
    snakemake.output.data, index=False, compression="gzip"
)
