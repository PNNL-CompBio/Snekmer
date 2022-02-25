"""search.smk: Module for supervised modeling from kmer vectors.
JEM: Currently this only creates features for a search which
     must be implemented in a standalone script.

author: @christinehc, @biodataganache
"""
# fmt: off
# snakemake config
# include: "kmerize.smk"
# localrules: all, generate, vectorize
# force snakemake v6.0+ (required for modules)
from snakemake.utils import min_version
min_version("6.0")


# load modules
module process_input:
    snakefile: "process_input.smk"
    config: config


module kmerize:
    snakefile: "kmerize.smk"
    config: config


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

# external libraries
import snekmer as skm
import numpy as np
#import pandas as pd
import matplotlib.pyplot as plt
from pandas import DataFrame, read_csv, read_json
from Bio import SeqIO
from sklearn.linear_model import LogisticRegressionCV, LogisticRegression
from sklearn.model_selection import train_test_split, StratifiedKFold, cross_val_score
from sklearn.preprocessing import LabelEncoder
from sklearn.ensemble import RandomForestClassifier
from sklearn.tree import DecisionTreeClassifier

# change matplotlib backend to non-interactive
plt.switch_backend("Agg")

# collect all fasta-like files, unzipped filenames, and basenames
input_files = glob(join("input", "*"))
zipped = [fa for fa in input_files if fa.endswith(".gz")]
unzipped = [
    fa.rstrip(".gz")
    for fa, ext in product(input_files, config["input"]["file_extensions"])
    if fa.rstrip(".gz").endswith(f".{ext}")
]

# map extensions to basename (basename.ext.gz -> {basename: ext})
uz_map = {
    skm.utils.split_file_ext(f)[0]: skm.utils.split_file_ext(f)[1] for f in zipped
}
fa_map = {
    skm.utils.split_file_ext(f)[0]: skm.utils.split_file_ext(f)[1] for f in unzipped
}
UZS = list(uz_map.keys())
FAS = list(fa_map.keys())

# parse any background files
bg_files = glob(join("input", "background", "*"))
if len(bg_files) > 0:
    bg_files = [skm.utils.split_file_ext(basename(f))[0] for f in bg_files]
NON_BGS, BGS = [f for f in FAS if f not in bg_files], bg_files

# terminate with error if invalid alphabet specified
skm.alphabet.check_valid(config["alphabet"])

# define output directory (helpful for multiple runs)
out_dir = skm.io.define_output_dir(config['alphabet'], config['k'],
                                   nested=config['output']['nested_dir'])

# define output files to be created by snekmer
rule all:
    input:
        expand(join("input", '{uz}'), uz=UZS),  # require unzipping
        expand(join(out_dir, "features", "{nb}", "{fa}.json.gz"), nb=NON_BGS, fa=FAS),  # correctly build features
        expand(join(out_dir, "model", "results", "{nb}.csv"), nb=NON_BGS)  # require model-building


# if any files are gzip zipped, unzip them
use rule unzip from process_input with:
    output:
        join("input", '{uz}')  # or join("input", "{uz}.{uzext}") ?

# JEM:
# Removing these two rules for the 'search' function since they
#          shouldn't be necessary

# read and process parameters from config
#use rule preprocess from process_input with:
#    input:
#        fasta=lambda wildcards: join("input", f"{wildcards.nb}.{fa_map[wildcards.nb]}")
#    output:
#        data=join(out_dir, "processed", "{nb}.json"),
#        desc=join(out_dir, "processed", "{nb}_description.csv")
#    log:
#        join(out_dir, "processed", "log", "{nb}.log")


# generate kmer features space from user params
#use rule generate from kmerize with:
#    input:
#        params=join(out_dir, "processed", "{nb}.json")
#    output:
#        labels=join(out_dir, "labels", "{nb}.txt")
#    log:
#        join(out_dir, "labels", "log", "{nb}.log")


# build kmer count vectors for each basis set
# JEM: vectorize_search doesn't depend on preprocessing
#      for inputs - takes directly from configfile
use rule vectorize_search from kmerize with:
    input:
        #kmers=join(out_dir, "labels", "{nb}.txt"),
        #params=join(out_dir, "processed", "{nb}.json"),
        fastas=unzipped
    log:
        join(out_dir, "features", "log", "{nb}.log")
    output:
        files=expand(join(out_dir, "features", "{{nb}}", "{fa}.json.gz"),
                     fa=FAS)

# JEM: all I've done so far is to change the name of the rule
rule search:
    input:
        vecfile=rules.vectorize_search.output.files
    output:
        results=join(out_dir, "model", "results", "{nb}.csv")
        #figs=directory(join(out_dir, "model", "figures", "{nb}"))
    run:
        new_basis = []
        with open(config["input"]["feature_set"], "r") as f:
            new_basis = skm.io.read_output_kmers(config["input"]["feature_set"])

        #new_vec_file = "snekmer-app-4/output/features/input/input.json.gz"
        #new_basis_file = "snekmer-app-4/output/labels/input.txt"

        new_vec_df = read_json(str(input.vecfile))
        #new_basis = skm.io.read_output_kmers(new_basis_file)
        new_vecs = skm.utils.to_feature_matrix(new_vec_df['vector'].values)


        #model_file = "snekmer/output/model/TS.pkl"
        #scorer_file = "snekmer/output/score/TS.pkl"

        # JEM: for now we'll specify the model we want
        #      to run in the config file
        model_file = config["input"]["model_file"]
        score_file = config["input"]["score_file"]

        # load model
        with open(model_file, 'rb') as mf:
            model = pickle.load(mf)

        #load scorer
        with open(score_file, 'rb') as sf:
            scorer = pickle.load(sf)

        new_scores = scorer.predict(new_vecs, new_basis)
        predictions = model.predict(new_scores.reshape(-1,1))

        new_vec_df['predicted'] = [True if p==1 else False for p in predictions]

        predictions_proba = model.predict_proba(new_scores.reshape(-1,1))
        new_vec_df['probabilities'] = [p[1] for p in predictions_proba]

        # save full results
        DataFrame(new_vec_df).to_csv(output.results, index=False)

# [in-progress] kmer walk
# if config['walk']:
    # use rule perform_kmer_walk from process_input with:
        # output:

# SUPERVISED WORKFLOW
# rule score:
#     input:
#         kmers=join(out_dir, "labels", "{nb}.txt"),
#         files=expand(join(out_dir, "features", "{{nb}}", "{fa}.json.gz"),
#                      fa=FAS)
#     output:
#         df=join(out_dir, "features", "score", "{nb}.csv.gz"),
#         scores=join(out_dir, "score", "weights", "{nb}.csv.gz"),
#         scorer=join(out_dir, "score", "{nb}.pkl")
#     log:
#         join(out_dir, "score", "log", "{nb}.log")
#     run:
#         # log script start time
#         start_time = datetime.now()
#         with open(log[0], "a") as f:
#             f.write(f"start time:\t{start_time}\n")
#
#         # get kmers for this particular set of sequences
#         kmers = skm.io.read_output_kmers(input.kmers)
#
#         # parse all data and label background files
#         label = config["score"]["lname"]
#         data = skm.io.vecfiles_to_df(
#             input.files, labels=config["score"]["labels"], label_name=label
#         )
#         data["background"] = [skm.utils.split_file_ext(f)[0] in BGS for f in data["filename"]]
#
#         # log conversion step runtime
#         skm.utils.log_runtime(log[0], start_time, step="vecfiles_to_df")
#
#         # parse family names and only add if some are valid
#         families = [
#             skm.utils.get_family(
#                 skm.utils.split_file_ext(fn)[0], regex=config["input"]["regex"]
#             )
#             for fn in data["filename"]
#         ]
#         if any(families):
#             label = "family"
#             data[label] = families
#
#         # binary T/F for classification into family
#         family = skm.utils.get_family(wildcards.nb)
#         binary_labels = [True if value == family else False for value in data[label]]
#
#         # define k-fold split indices
#         if config["model"]["cv"] > 1:
#             cv = StratifiedKFold(n_splits=config["model"]["cv"], shuffle=True)
#
#             # stratify splits by [0,1] family assignment
#             for n, (i_train, _) in enumerate(cv.split(data["vector"], binary_labels)):
#                 data[f"train_cv-{n + 1:02d}"] = [idx in i_train for idx in data.index]
#
#         elif config["model"]["cv"] in [0, 1]:
#             i_train, _ = train_test_split(data.index, stratify=binary_labels)
#             data["train"] = [idx in i_train for idx in data.index]
#
#         # generate family scores and object
#         scorer = skm.model.KmerScorer()
#         scorer.fit(
#             kmers,
#             data,
#             skm.utils.get_family(wildcards.nb, regex=config["input"]["regex"]),
#             label_col=label,
#             **config["score"]["scaler_kwargs"],
#         )
#
#         # append scored sequences to dataframe
#         data = data.merge(DataFrame(scorer.scores["sample"]), left_index=True, right_index=True)
#         if data.empty:
#             raise ValueError("Blank df")
#
#         # save score loadings
#         class_probabilities = (
#             DataFrame(scorer.probabilities, index=scorer.kmers.basis)
#             .reset_index()
#             .rename(columns={"index": "kmer"})
#         )
#
#         # log time to compute class probabilities
#         skm.utils.log_runtime(log[0], start_time, step="class_probabilities")
#
#         # save all files to respective outputs
#         delete_cols = ["vec", "vector"]
#         for col in delete_cols:
#             # if col in data.columns:
#             #     data = data.drop(columns=col)
#             if col in class_probabilities.columns:
#                 class_probabilities = class_probabilities.drop(columns=col)
#         data.to_csv(output.df, index=False, compression="gzip")
#         class_probabilities.to_csv(output.scores, index=False, compression="gzip")
#         with open(output.scorer, "wb") as f:
#             pickle.dump(scorer, f)
#
#         # record script endtime
#         skm.utils.log_runtime(log[0], start_time)
#
