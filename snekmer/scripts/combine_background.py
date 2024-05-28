# ---------------------------------------------------------
# Imports
# ---------------------------------------------------------
from datetime import datetime

import numpy as np
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
    # load input files
    kmerlist, df = skm.io.load_npz(bgf)
    kmerlist = np.hstack(kmerlist)

    # create matrix of kmer counts
    counts = np.vstack(df["sequence_vector"].values)
    x, y = len(df), len(counts[0])
    del df  # clear memory

    # sum counts per kmer and harmonize
    counts = np.sum(counts, axis=0)
    counts = basis.harmonize(counts, kmerlist)

    # combine counts with full matrix
    total_norm += x
    matrix = np.add(matrix, counts)

    del counts

# normalize by total # background seqs
matrix = np.divide(matrix, total_norm)

# save new background kmer vector
np.savez_compressed(
    snakemake.output[0],
    kmer_list=basis.basis.basis,
    kmer_counts=matrix,
)
