# force snakemake v6.0+ (required for modules)
from snakemake.utils import min_version

min_version("6.0")


# load snakemake modules
module process_input:
    snakefile:
        "process_input.smk"
    config:
        config


module kmerize:
    snakefile:
        "kmerize.smk"
    config:
        config


# built-in imports
import gzip
import json
import pickle
from datetime import datetime
from glob import glob
from itertools import product, repeat
from multiprocessing import Pool
from os import makedirs
from os.path import basename, dirname, exists, join, splitext

import matplotlib.pyplot as plt
import numpy as np
from Bio import SeqIO
from pandas import DataFrame, read_csv, read_json
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegressionCV
from sklearn.model_selection import (StratifiedKFold, cross_val_score,
                                     train_test_split)
from sklearn.preprocessing import LabelEncoder
from sklearn.tree import DecisionTreeClassifier

# external libraries
import snekmer as skm

# change matplotlib backend to non-interactive
plt.switch_backend("Agg")

# collect all fasta-like files, unzipped filenames, and basenames
input_dir = "input" if str(config["input_dir"]) == "None" else config["input_dir"]
input_files = glob(join(input_dir, "*"))
zipped = [fa for fa in input_files if fa.endswith(".gz")]
unzipped = [
    fa.rstrip(".gz")
    for fa, ext in product(input_files, config["input"]["file_extensions"])
    if fa.rstrip(".gz").endswith(f".{ext}")
]

# map extensions to basename (basename.ext.gz -> {basename: ext})
UZ_MAP = {
    skm.utils.split_file_ext(f)[0]: skm.utils.split_file_ext(f)[1] for f in zipped
}
FA_MAP = {
    skm.utils.split_file_ext(f)[0]: skm.utils.split_file_ext(f)[1] for f in unzipped
}

# get unzipped filenames
UZS = [f"{f}.{ext}" for f, ext in UZ_MAP.items()]

# isolate basenames for all files
FAS = list(FA_MAP.keys())

# parse any background files
bg_files = glob(join(input_dir, "background", "*"))
if len(bg_files) > 0:
    bg_files = [skm.utils.split_file_ext(basename(f))[0] for f in bg_files]
NON_BGS, BGS = [f for f in FAS if f not in bg_files], bg_files

# terminate with error if invalid alphabet specified
skm.alphabet.check_valid(config["alphabet"])

# define output directory (helpful for multiple runs)
out_dir = skm.io.define_output_dir(
    config["alphabet"], config["k"], nested=config["output"]["nested_dir"]
)


# define output files to be created by snekmer
rule all:
    input:
        expand(join(input_dir, "{uz}"), uz=UZS),  # require unzipping
        expand(join(out_dir, "cluster", "{nb}.pkl"), nb=NON_BGS),  # require model-building


# if any files are gzip zipped, unzip them
use rule unzip from process_input with:
    output:
        join(input_dir, "{uz}"),


# read and process parameters from config
use rule preprocess from process_input with:
    input:
        fasta=lambda wildcards: join(input_dir, f"{wildcards.nb}.{FA_MAP[wildcards.nb]}"),
    output:
        data=join(out_dir, "processed", "full", "{nb}.json"),
        desc=join(out_dir, "processed", "full", "{nb}_description.csv"),
    log:
        join(out_dir, "processed", "log", "{nb}.log")


# generate kmer features space from user params
# use rule generate from kmerize with:
#     input:
#         params=join(out_dir, "processed", "full", "{nb}.json"),
#     output:
#         labels=join(out_dir, "labels", "full", "{nb}.txt"),
#     log:
#         join(out_dir, "labels", "log", "{nb}.log")


# build kmer count vectors for each basis set
use rule vectorize_full from kmerize with:
    input:
        params=join(out_dir, "processed", "{nb}.json"),
        fasta=lambda wildcards: join(input_dir, f"{wildcards.nb}.{FA_MAP[wildcards.nb]}")
    log:
        join(out_dir, "features", "log", "{nb}.log")
    output:
        file=join(out_dir, "features", "full", "{nb}.json.gz"),
        kmers=join(out_dir, "labels", "full", "{nb}.txt")


# [in-progress] kmer walk
# if config['walk']:
# use rule perform_kmer_walk from process_input with:
# output:


# UNSUPERVISED WORKFLOW
rule cluster:
    input:
        files=expand(join(out_dir, "features", "full", "{fa}.json.gz"), fa=NON_BGS)
    output:
        model=join(out_dir, "cluster", "{nb}.pkl"),
        figs=directory(join(out_dir, "cluster", "figures", "{nb}")),
    log:
        join(out_dir, "cluster", "log", "{nb}.log"),
    run:
        # log script start time
        start_time = datetime.now()
        with open(log[0], "a") as f:
            f.write(f"start time:\t{start_time}\n")

        # parse all data and label background files
        label = config["score"]["lname"]
        data = skm.io.vecfiles_to_df(
            input.files, labels=config["score"]["labels"], label_name=label
        )
        data["background"] = [
            skm.utils.split_file_ext(f)[0] in BGS for f in data["filename"]
        ]

        # log conversion step runtime
        skm.utils.log_runtime(log[0], start_time, step="vecfiles_to_df")

        # define feature matrix of kmer vectors not from background set
        bg, non_bg = data[data["background"]], data[~data["background"]]
        full_feature_matrix = skm.utils.to_feature_matrix(data["vector"].values)
        feature_matrix = skm.utils.to_feature_matrix(non_bg["vector"].values)

        # currently not used
        if len(bg) > 0:
            bg_feature_matrix = skm.utils.to_feature_matrix(bg["vector"].values)

        # fit and save clustering model
        model = skm.cluster.KmerClustering(
            config["cluster"]["method"], config["cluster"]["params"]
        )
        model.fit(full_feature_matrix)
        with open(output.model, "wb") as f:
            pickle.dump(model, f)

        # log time to compute clusters
        skm.utils.log_runtime(log[0], start_time, step="clustering")

        # insert plots here?
        if not exists(output.figs):
            makedirs(output.figs)
        fig, ax = skm.plot.show_explained_variance_curve(full_feature_matrix)
        fig.savefig(join(output.figs, "pca_explained_variance_curve.png"))
        plt.close("all")

        fig, ax = skm.plot.get_tsne_clusters(full_feature_matrix, model.model.labels_)
        fig.savefig(join(output.figs, "tsne_clusters.png"))
        plt.close("all")

        # record script endtime
        skm.utils.log_runtime(log[0], start_time)
