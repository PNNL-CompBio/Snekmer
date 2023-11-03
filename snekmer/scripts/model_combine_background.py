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

# load all background files to build combined bases
kmers, scores, norms = list(), list(), list(), list()
for f in snakemake.input.scores:
    kmerlist_bg, probas_bg, scores_bg = skm.io.load_npz(bgf)
    kmers.append(kmerlist_bg)
    norms.append(round(1 / probas_bg["probability"].min()))
    scores.append(scores_bg)

# rescale probabilities to respective file sizes (nseqs)
scaled = list()
for i, p in enumerate(scores):
    scaled.append(s * norms[i])

# make a superset of kmers
basis = skm.vectorize.KmerVec(config["alphabet"], config["k"])
basis.set_kmer_set(np.unique(np.concatenate(kmers)))

# rescale and harmonize vecs before combining
score_vec = np.zeros(len(np.hstack(kmers)))
for i, s in enumerate(scaled):
    converted = basis.harmonize(s.reshape(-1, 1).T, kmers[i])
    score_vec += converted

# renormalize combined scores from all bg files w/ total nseqs
total_norm = sum(norms)
score_vec = score_vec / total_norm

# save new background score vector
np.savez_compressed(
    output[0],
    kmerlist=kmers,
    scores=score_vec,
)
