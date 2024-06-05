"""
Created on Tue Apr 25 15:14:39 2023

author: @tnitka
"""

# ---------------------------------------------------------
# Imports
# ---------------------------------------------------------
import pickle

import snekmer as skm
import pandas as pd
import numpy as np
import gzip
import gc

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

with gzip.open(snakemake.input.weights, "rb") as f:
    weights = pd.read_csv(f)

with gzip.open(snakemake.input.scores, "rb") as f:
    scores = pd.read_csv(f)
    scores = scores.astype("float64")
    scores = scores.to_numpy()


with gzip.open(snakemake.input.kmers, "rb") as f:
    kmers = pd.read_csv(f)

# set category label name (e.g. "family")
label = config["score"]["lname"] if str(config["score"]["lname"]) != "None" else "label"

if config["k"] == 2:
    kmers = kmers.fillna("NA")
    scores = scores.fillna("NA")

family = skm.utils.get_family(
    skm.utils.split_file_ext(snakemake.input.weights)[0],
    regex=config["input_file_regex"],
)
scorer = skm.score.KmerScorer()

gc.collect()

# set number of permutations to test
n_iter = config["motif"]["n"]


# binary T/F for classification into family
family = skm.utils.get_family(snakemake.wildcards.nb)
binary_labels = [True if value == family else False for value in data[label]]

# get alphabet name
if config["alphabet"] in skm.alphabet.ALPHABET_ORDER.keys():
    alphabet_name = skm.alphabet.ALPHABET_ORDER[config["alphabet"]].capitalize()
else:
    alphabet_name = str(config["alphabet"]).capitalize()


# gather scores from rescoring jobs
score_matrix = kmers.rename(columns={"0": "kmer"})
score_array = pd.DataFrame.to_numpy(score_matrix)
motif = skm.motif.SnekmerMotif()
for file in snakemake.input.perm_scores:
    with gzip.open(file) as f:
        perm_scores = pd.DataFrame.to_numpy(pd.read_csv(f))
    score_array = np.hstack((score_array, perm_scores))

score_array = np.delete(score_array, 0, 1)
score_matrix = score_matrix.merge(
    pd.DataFrame(score_array), left_index=True, right_index=True
)

# format output
scores = np.ravel(scores)
output_matrix = motif.p_values(score_matrix, scores, n_iter)
output_matrix = output_matrix.astype(
    {
        "kmer": "str",
        "real score": "float32",
        "false positives": "int32",
        "n": "int32",
        "p": "float32",
    }
)
output_matrix.sort_values(by=["p", "real score"], ascending=[True, False], inplace=True)

# save output
score_matrix.to_csv(snakemake.output.data, index=False, compression="gzip")
output_matrix.to_csv(snakemake.output.p_values, index=False)
