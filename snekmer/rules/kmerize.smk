"""kmerize.smk: Module for kmer vector generation.

author: @christinehc

"""


# include unzipping module
include: "process.smk"


# built-in imports
import itertools
import gzip
import json
import pickle
from datetime import datetime
from glob import glob
from os.path import basename, exists, join

# external libraries
import numpy as np
import pandas as pd
import snekmer as skm
from Bio import SeqIO

# get input files
input_dir = (
    "input"
    if (("input_dir" not in config) or (str(config["input_dir"]) == "None"))
    else config["input_dir"]
)
input_files = glob(join(input_dir, "*"))

input_file_exts = ["fasta", "fna", "faa", "fa"]
if "input_file_exts" in config:
    input_file_exts = config["input_file_exts"]

unzipped = [
    fa.rstrip(".gz")
    for fa, ext in itertools.product(input_files, input_file_exts)
    if fa.rstrip(".gz").endswith(f".{ext}")
]
zipped = [fa for fa in input_files if fa.endswith(".gz")]
UZS = [skm.utils.split_file_ext(f)[0] for f in zipped]
FAS = [skm.utils.split_file_ext(f)[0] for f in unzipped]

# map extensions to basename (basename.ext.gz -> {basename: ext})
UZ_MAP = {
    skm.utils.split_file_ext(f)[0]: skm.utils.split_file_ext(f)[1] for f in zipped
}
FA_MAP = {
    skm.utils.split_file_ext(f)[0]: skm.utils.split_file_ext(f)[1] for f in unzipped
}


rule vectorize:
    input:
        fasta=lambda wildcards: join(
            "input", f"{wildcards.nb}.{FA_MAP[wildcards.nb]}"
        ),
        kmerbasis=join(input_dir, "basis.txt"), # this is optional
    output:
        data=join("output", "vector", "{nb}.npz"),
        kmerobj=join("output", "kmerize", "{nb}.kmers"),
    log:
        join("output", "kmerize", "log", "{nb}.log"),
    run:
        # initialize kmerization object
        alphabet = config["alphabet"]
        k = config["k"]

        kmer = skm.vectorize.KmerVec(alphabet=alphabet, k=k)
        # we want to check to see if this is a ridiculous combination

        # read kmerbasis if present
        min_filter = 0
        if hasattr(input, "kmerbasis") and exists(input.kmerbasis):
            kmerbasis = skm.io.read_kmers(input.kmerbasis)

            # quick way to get the number of proteins in
            #     the fasta file so we can set up an array
            #     ahead of time
            nprot = len([1 for line in open(input.fasta) if line.startswith(">")])

        else:
            # we make our own kmerbasis and filter for minimum
            #    number of occurrences, etc.
            # we will only allow filtering by number of kmers
            # if we're not using a basis set as input -
            if "min_filter" in config:
                min_filter = config["min_filter"]

            # make basis
            kmerbasis = {}
            fasta = SeqIO.parse(input.fasta, "fasta")

            nprot = 0
            for f in fasta:
                nprot += 1
                these = kmer.reduce_vectorize(f.seq)
                for key in these:
                    if key in kmerbasis:
                        kmerbasis[key] += 1
                    else:
                        kmerbasis[key] = 1

            # if min_filter is specified as less than 1 it
            #    will be treated as a percentage of the total
            #    number of kmers *present* in the input set
            # And I think it makes sense to have this be in the
            #    reverse sense as the integer min_filter. That is
            #    we will filter out 1-min_filter percentage of kmers
            if min_filter > 0 and min_filter < 1:
                min_filter = len(list(kmerbasis.keys()))*(1.0-min_filter)

            kmerbasis = np.array(list(kmerbasis.keys()))[
                np.array(list(kmerbasis.values())) > min_filter
            ]

            # an issue here is that the filter might remove all kmers
            #    from the basis and break the code. We should check for
            # this and exit gracefully
            # So for now we will die and write a helpful message to
            #  the logfile
            if len(kmerbasis) == 0:
                msg = "FATAL: min_filter application results in no remaining kmers\n"
                print(msg)
                with open(log[0], "a") as f:
                    f.write(msg)
                # this state will cause Snekmer to crash downstream

        kmer.set_kmer_set(kmerbasis)

        # (re)read fasta using bioconda obj
        fasta = SeqIO.parse(input.fasta, "fasta")

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
                    alphabet=config["alphabet"],
                    mapping=skm.alphabet.FULL_ALPHABETS,
                )
            )
            ids.append(f.id)
            lengths.append(len(f.seq))

        # save seqIO output and transformed vecs
        np.savez_compressed(
            output.data,
            kmerlist=kmerbasis,
            ids=ids,
            seqs=seqs,
            vecs=vecs,
            lengths=lengths,
        )

        with open(output.kmerobj, "wb") as f:
            pickle.dump(kmer, f)
