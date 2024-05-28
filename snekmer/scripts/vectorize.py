# imports
import pickle
from os.path import exists

import numpy as np
import snekmer as skm
from Bio import SeqIO

# run script
kmer = skm.vectorize.KmerVec(
    alphabet=snakemake.config["alphabet"], k=snakemake.config["k"]
)

# read kmerbasis if present
min_filter = 0
if hasattr(snakemake.input, "kmerbasis") and exists(snakemake.input.kmerbasis):
    kmerbasis = skm.io.read_kmers(snakemake.input.kmerbasis)

    # quick way to get the number of proteins in
    #     the fasta file so we can set up an array
    #     ahead of time
    nprot = len([1 for line in open(snakemake.input.fasta) if line.startswith(">")])

else:
    # we make our own kmerbasis and filter for minimum
    #    number of occurrences, etc.
    # we will only allow filtering by number of kmers
    # if we're not using a basis set as input -
    if "min_filter" in snakemake.config:
        min_filter = snakemake.config["min_filter"]

        # make basis
    kmerbasis = {}
    fasta = SeqIO.parse(snakemake.input.fasta, "fasta")

    nprot = 0
    for f in fasta:
        nprot += 1
        these = kmer.reduce_vectorize(f.seq)
        for key in these:
            if key in kmerbasis:
                kmerbasis[key] += 1
            else:
                kmerbasis[key] = 1

    kmerbasis = np.array(list(kmerbasis.keys()))[
        np.array(list(kmerbasis.values())) > min_filter
    ]

kmer.set_kmer_set(kmerbasis)

# (re)read fasta using bioconda obj
fasta = SeqIO.parse(snakemake.input.fasta, "fasta")

# pre-allocate an array to keep results
vecs = np.zeros((nprot, len(kmerbasis)))

# I question whether we need to keep the reduced seqs here
seqs, ids, lengths = list(), list(), list()
n = 0
for f in fasta:
    addvec = kmer.reduce_vectorize(f.seq)
    vecs[n][np.isin(kmerbasis, addvec)] = 1
    n += 1
    seqs.append(
        skm.vectorize.reduce(
            f.seq,
            alphabet=snakemake.config["alphabet"],
            mapping=skm.alphabet.FULL_ALPHABETS,
        )
    )
    ids.append(f.id)
    lengths.append(len(f.seq))

    # save seqIO output and transformed vecs
np.savez_compressed(
    snakemake.output.data,
    kmerlist=kmerbasis,
    ids=ids,
    seqs=seqs,
    vecs=vecs,
    lengths=lengths,
)

with open(snakemake.output.kmerobj, "wb") as f:
    pickle.dump(kmer, f)
