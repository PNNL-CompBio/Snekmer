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

# load kmer basis for family of interest
basis = skm.io.load_pickle(snakemake.input.kmerobj)

# load all background files, rescale, and harmonize vecs
score_vec = np.zeros(len(np.hstack(basis.basis.basis)))
total_norm = 0
for bgf in snakemake.input.scores:
    loaded = np.load(bgf, allow_pickle=True)
    weight = config["score"]["background_weight"]
    scores = weight * loaded["probas"]

    norm = round(1 / loaded["probas"].min())
    scaled = scores * norm
    harmonized = basis.harmonize(scaled.reshape(-1, 1).T, loaded["kmerlist"][0])

    total_norm += norm
    score_vec = np.sum([score_vec, harmonized], axis=0)

# rescale probabilities to respective file sizes (nseqs)
# scaled = list()
# for i, p in enumerate(scores):
#     scaled.append(s * norms[i])

# make a superset of kmers
# basis = skm.vectorize.KmerVec(config["alphabet"], config["k"])
# basis.set_kmer_set(np.unique(np.concatenate(kmers)))

# rescale vecs before combining
# score_vec = np.zeros(len(np.hstack(basis.basis.basis)))
# for i, s in enumerate(scaled):
#     # converted = basis.harmonize(s.reshape(-1, 1).T, kmers[i])
#     score_vec += s

# renormalize combined scores from all bg files w/ total nseqs
# total_norm = sum(norms)
score_vec = score_vec / total_norm

# save new background score vector
np.savez_compressed(
    snakemake.output[0],
    kmerlist=basis.basis.basis,
    scores=score_vec,
)
