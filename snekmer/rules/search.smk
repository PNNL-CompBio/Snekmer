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
        expand(join("input", "{uz}"), uz=UZS),  # require unzipping
        expand(join(out_dir, "vector", "{f}.npz"), f=FILES),
        expand(join(out_dir, "search", "{fam}.csv"), fam=FAMILIES),  # require search


# if any files are gzip zipped, unzip them
if len(UZS) > 0:
    use rule unzip from process with:
        output:
            join("input", "{uz}"),  # or join("input", "{uz}.{uzext}") ?
# build kmer count vectors for each basis set


use rule vectorize from kmerize with:
    input:
        fasta=lambda wildcards: join("input", f"{wildcards.f}.{FILE_MAP[wildcards.f]}"),
    output:
        data=join(out_dir, "vector", "{f}.npz"),
        kmerobj=join(out_dir, "kmerize", "{f}.kmers"),
    log:
        join(out_dir, "kmerize", "log", "{f}.log"),


rule search:
    input:
        files=expand(join(out_dir, "vector", "{f}.npz"), f=FILES),  # change to data=join("output", "vector", "{nb}.npz")
        model=join(config["model_dir"], "{fam}.model"),
        kmerobj=join(config["basis_dir"], "{fam}.kmers"),
        scorer=join(config["score_dir"], "{fam}.scorer"),
    output:
        results=join(out_dir, "search", "{fam}.csv"),
    run:
        # simplify variable name
        family = wildcards.fam

        # get kmers for this particular set of sequences
        kmer = skm.io.load_pickle(input.kmerobj)
        model = skm.io.load_pickle(input.model)
        scorer = skm.io.load_pickle(input.scorer)

        # load vectorized sequences, score, and predict scores
        results = list()
        for f in input.files:
            df, kmerlist = skm.io.load_npz(f)
            filename = skm.utils.split_file_ext(basename(f))[0]

            vecs = skm.utils.to_feature_matrix(df["sequence_vector"].values)
            #vecs = df["sequence_vector"].values
            #print(vecs)
            base_kmers = list(kmer.kmer_set.kmers)
            scores = scorer.predict(vecs, kmerlist)

            predictions = model.predict(scores.reshape(-1, 1))
            predicted_probas = model.predict_proba(scores.reshape(-1, 1))
            
            # display results (score, family assignment, and probability)
            df[f"{family}_score"] = scores  # scorer output
            df[family] = [True if p == 1 else False for p in predictions]
            df[f"{family}_probability"] = [p[1] for p in predicted_probas]
            df["filename"] = f"{filename}.{FILE_MAP[filename]}"
            #df["seq_length"] = seq_length
            #df["seq_ids"] = seq_id

            results.append(df)

        # save full results
        results = pd.concat(results, ignore_index=True).drop(columns=["sequence_vector"])
        results["model"] = basename(input.model)
        results.to_csv(output.results, index=False)
