# ---------------------------------------------------------
# Imports
# ---------------------------------------------------------
from datetime import datetime

import numpy as np
import pandas as pd
import snekmer as skm

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


# # get kmers for this particular set of sequences
# kmers = skm.io.load_pickle(snakemake.input.kmerobj)

# # load vectorized seq data
# kmerlist, df = skm.io.load_npz(snakemake.input.data)

# # create matrix of kmer counts
# x, y = len(df["sequence_vector"]), len(df["sequence_vector"][0])
# matrix = np.zeros(x * y).reshape((x, y))
# for i in range(x):
#     for j in range(y):
#         value = df["sequence_vector"][i]
#         value = value[j]
#         matrix[i, j] = value

# probas = skm.score.feature_class_probabilities(
#     matrix.T, df["filename"], kmers=kmers.basis.basis
# )

# # save score weights
# np.savez_compressed(
#     snakemake.output.scores,
#     kmerlist=kmerlist,
#     probas=probas["probability"].values,
# )


# load all background files and harmonize vecs
matrix = np.zeros(len(np.hstack(basis.basis.basis)))
total_norm = 0
for bgf in snakemake.input.data:
    kmerlist, df = skm.io.load_npz(bgf)

    # create matrix of kmer counts
    x, y = len(df), len(df["sequence_vector"].to_numpy()[0])
    counts = np.zeros(x * y).reshape((x, y))
    for i in range(x):
        for j in range(y):
            value = df["sequence_vector"].to_numpy()[i]
            value = value[j]
            counts[i, j] = value

    harmonized = basis.harmonize(counts, np.hstack(kmerlist))
    summed = np.sum(harmonized, axis=0)

    total_norm += x
    matrix = np.add(matrix, summed)

# normalize by total # background seqs
matrix = matrix / total_norm

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
# score_vec = score_vec / total_norm

# save new background kmer vector
np.savez_compressed(
    snakemake.output[0],
    kmer_list=basis.basis.basis,
    kmer_counts=matrix,
)
