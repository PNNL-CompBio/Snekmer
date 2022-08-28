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
from os.path import basename, join

# external libraries
import numpy as np
import snekmer as skm
from Bio import SeqIO

# get input files
input_dir = "input" if (("input_dir" not in config) or (str(config["input_dir"]) == "None")) else config["input_dir"]
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
    output:
        data=join("output", "vector", "{nb}.npz"),
        kmerobj=join("output", "kmerize", "{nb}.kmers"),
    log:
        join("output", "kmerize", "log", "{nb}.log"),
    run:
        # read fasta using bioconda obj
        fasta = SeqIO.parse(input.fasta, "fasta")

        # initialize kmerization object
        kmer = skm.vectorize.KmerVec(alphabet=config["alphabet"], k=config["k"])

        vecs, seqs, ids, lengths = list(), list(), list(), list()
        for f in fasta:
            vecs.append(kmer.reduce_vectorize(f.seq))
            seqs.append(
                skm.vectorize.reduce(
                    f.seq,
                    alphabet=config["alphabet"],
                    mapping=skm.alphabet.FULL_ALPHABETS,
                )
            )
            ids.append(f.id)
            lengths.append(len(f.seq))

        # memfix: we need to reformat the vec array to reflect
        #         the all observed kmers in a matrix - rather than
        #         a list of lists of variable length. Then we'll also
        #         add the kmer order to kmer to what we save - that way we
        #         don't need to consider every possible kmer
        vecs, kmerlist = skm.vectorize.make_feature_matrix(vecs)

        kmer.set_kmer_set(kmer_set=kmerlist)

        # save seqIO output and transformed vecs
        np.savez_compressed(output.data, kmerlist=kmerlist, ids=ids,
                        seqs=seqs, vecs=vecs, lengths=lengths)

        with open(output.kmerobj, "wb") as f:
            pickle.dump(kmer, f)
