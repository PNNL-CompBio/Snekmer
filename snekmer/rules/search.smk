"""search.smk: Module for supervised modeling from kmer vectors.

author: @christinehc, @biodataganache

"""
# snakemake config
from snakemake.utils import min_version

min_version("6.0")  # force snakemake v6.0+ (required for modules)
ruleorder: vectorize > search


# load modules
module process:
    snakefile:
        "process.smk"
    config:
        config

module kmerize:
    snakefile:
        # "snekmer.smk"
        "kmerize.smk"
    config:
        config


# built-in imports
import gzip
import json
import pickle
from ast import literal_eval
from datetime import datetime
from glob import glob
from itertools import product, repeat
from multiprocessing import Pool
from os import makedirs
from os.path import basename, dirname, exists, join, splitext

# import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import SeqIO
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression, LogisticRegressionCV
from sklearn.model_selection import (StratifiedKFold, cross_val_score,
                                     train_test_split)
from sklearn.preprocessing import LabelEncoder
from sklearn.tree import DecisionTreeClassifier

# external libraries
import snekmer as skm

# change matplotlib backend to non-interactive
plt.switch_backend("Agg")

# collect all fasta-like files, unzipped filenames, and basenames
input_dir = "input" if (("input_dir" not in config) or (str(config["input_dir"]) == "None")) else config["input_dir"]
input_files = glob(join(input_dir, "*"))

model_files = glob(join(config["model_dir"], "*.model"))
zipped = [f for f in input_files if f.endswith(".gz")]

input_file_exts = ["fasta", "fna", "faa", "fa"]
if "input_file_exts" in config:
    input_file_exts = config["input_file_exts"]

unzipped = [
    f.rstrip(".gz")
    for f, ext in product(input_files, input_file_exts)
    if f.rstrip(".gz").endswith(f".{ext}")
]

# map extensions to basename (basename.ext.gz -> {basename: ext})
UZ_MAP = {
    skm.utils.split_file_ext(f)[0]: skm.utils.split_file_ext(f)[1] for f in zipped
}
FILE_MAP = {
    skm.utils.split_file_ext(f)[0]: skm.utils.split_file_ext(f)[1] for f in unzipped
}
UZS = [f"{f}.{ext}" for f, ext in UZ_MAP.items()]
FILES = list(FILE_MAP.keys())

input_file_regex = ".*"

FAMILIES = [skm.utils.get_family(f, regex=config["input_file_regex"]) for f in model_files]

# define output directory (helpful for multiple runs)
out_dir = skm.io.define_output_dir(
    config["alphabet"], config["k"], nested=config["nested_output"]
)

# define output files to be created by snekmer
rule all:
    input:
        expand(join(input_dir, "{uz}"), uz=UZS),  # require unzipping
        join(config["basis_dir"], "search_kmers.txt"), # require common basis
        expand(join(out_dir, "vector", "{f}.npz"), f=FILES),
        expand(join(out_dir, "search", "{fam}", "{f}.csv"), fam=FAMILIES, f=FILES),   # require search


# if any files are gzip zipped, unzip them
if len(UZS) > 0:
    use rule unzip from process with:
        output:
            join("input", "{uz}"),  # or join("input", "{uz}.{uzext}") ?

# build kmer count vectors for each basis set
rule common_basis:
    input:
        kmerobjs=expand(join(config["basis_dir"], "{fam}.kmers"), fam=FAMILIES),
    output:
        kmerbasis=join(config["basis_dir"], "search_kmers.txt"),
    log:
        join(config["basis_dir"], "log", "common_basis.log"),
    run:
        common_basis = list()
        for kobj in input.kmerobjs:
            kmers = skm.io.load_pickle(kobj)

            # consolidate overly long lists of duplicate kmers
            if len(common_basis) > 1e10:
                common_basis = list(set(common_basis))
            common_basis.extend(list(kmers.kmer_set.kmers))

        # capture common basis set -- is faster than np.unique
        common_basis = set(common_basis)
        common_basis = sorted(list(common_basis))

        # output type: plaintext (no pandas) would likely be more compact
        with open(output.kmerbasis, "w") as f:
            for kmer in common_basis:
                f.write(f"{kmer}\n")
        # df = pd.DataFrame({'common': common_basis})
        # df.to_csv(output.kmerbasis, index=False)


use rule vectorize from kmerize with:
    input:
        fasta=lambda wildcards: join(input_dir, f"{wildcards.f}.{FILE_MAP[wildcards.f]}"),
        kmerbasis=rules.common_basis.output.kmerbasis,
    output:
        data=join(out_dir, "vector", "{f}.npz"),
        kmerobj=join(out_dir, "kmerize", "{f}.kmers"),
    log:
        join(out_dir, "kmerize", "log", "{f}.log"),


rule search:
    input:
        vecs=join(out_dir, "vector", "{f}.npz"),  # change to data=join("output", "vector", "{nb}.npz")
        model=join(config["model_dir"], "{fam}.model"),
        kmerobj=join(config["basis_dir"], "{fam}.kmers"),
        scorer=join(config["score_dir"], "{fam}.scorer"),
    output:
        results=join(out_dir, "search", "{fam}", "{f}.csv"),
    run:
        # simplify variable name
        family = wildcards.fam

        # get kmers for this particular set of sequences
        # print(f"starting {family}")
        kmer = skm.io.load_pickle(input.kmerobj)
        model = skm.io.load_pickle(input.model)
        scorer = skm.io.load_pickle(input.scorer)
        # print(f"loaded model {family}")

        # load vectorized sequences, score, and predict scores
        kmerlist, df  = skm.io.load_npz(input.vecs)
        filename = skm.utils.split_file_ext(basename(input.vecs))[0]

        # print(f"making feature matrix {family}")
        vecs = skm.utils.to_feature_matrix(df["sequence_vector"].values)

        # print(f"getting scores {family}")
        scores = scorer.predict(vecs, kmerlist[0])
        # print(f"making predictions {family}")
        predictions = model.predict(scores.reshape(-1, 1))
        # print(f"getting probabilities {family}")
        predicted_probas = model.predict_proba(scores.reshape(-1, 1))

        # display results (score, family assignment, and probability)
        df["score"] = scores  # scorer output
        df["in_family"] = [True if p == 1 else False for p in predictions]
        df["probability"] = [p[1] for p in predicted_probas]
        df["filename"] = f"{filename}.{FILE_MAP[filename]}"
        df["model"] = basename(input.model)

        df = df.drop(columns=["sequence_vector"])
        df.to_csv(output.results, index=False)
