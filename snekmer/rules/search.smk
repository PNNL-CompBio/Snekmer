"""search.smk: Module for supervised modeling from kmer vectors.
JEM: Currently this only creates features for a search which
     must be implemented in a standalone script.

author: @christinehc, @biodataganache
"""
# snakemake config
from snakemake.utils import min_version
min_version("6.0") # force snakemake v6.0+ (required for modules)


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
# import pandas as pd
import matplotlib.pyplot as plt
import pandas as pd

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
model_files = glob(join(config["model_dir"], "*.pkl"))
zipped = [f for f in input_files if f.endswith(".gz")]
unzipped = [
    f.rstrip(".gz")
    for f, ext in product(input_files, config["file_extensions"])
    if f.rstrip(".gz").endswith(f".{ext}")
]

# map extensions to basename (basename.ext.gz -> {basename: ext})
uz_map = {
    skm.utils.split_file_ext(f)[0]: skm.utils.split_file_ext(f)[1] for f in zipped
}
FILE_MAP = {
    skm.utils.split_file_ext(f)[0]: skm.utils.split_file_ext(f)[1] for f in unzipped
}
UZS = [f"{f}.{ext}" for f, ext in uz_map.items()]
FILES = list(FILE_MAP.keys())
FAMILIES = [skm.utils.get_family(f, regex=config["regex"]) for f in model_files]


# define output files to be created by snekmer
rule all:
    input:
        expand(join("input", "{uz}"), uz=UZS),  # require unzipping
        expand(join("output", "features", "{fam}", "{f}.json.gz"), fam=FAMILIES, f=FILES),  # correctly build features
        expand(join("output", "search", "{fam}.csv"), fam=FAMILIES)  # require model-building


# if any files are gzip zipped, unzip them
if len(UZS) > 0:
    use rule unzip from process_input with:
        output:
            join("input", "{uz}")  # or join("input", "{uz}.{uzext}") ?


# build kmer count vectors for each basis set
rule vectorize:
    input:
        fastas=unzipped,
        basis=join(config["basis_dir"], "{fam}.txt")
    log:
        join("output", "features", "log", "{fam}.log"),
    output:
        files=expand(
            join("output", "features", "{{fam}}", "{f}.json.gz"), f=FILES
        )
    run:
        start_time = datetime.now()

        # get kmers for this search
        kmers = skm.io.read_output_kmers(input.basis)

        # sort i/o lists to match wildcard order
        fastas = sorted(input.fastas)
        outfiles = sorted(output.files)

        # revectorize based on full kmer list
        for i, fa in enumerate(fastas):
            results = {"seq_id": [], "vector": []}
            seq_list, id_list = skm.io.read_fasta(fa)
            for seq, sid in zip(seq_list, id_list):
                results["seq_id"] += [sid]
                results["vector"] += [
                    skm.transform.vectorize_string(
                        seq,
                        config["k"],
                        config["alphabet"],
                        start=config["start"],
                        end=config["end"],
                        filter_list=kmers,
                        verbose=config["verbose"],
                        log_file=log[0],
                    )
                ]

            with gzip.open(outfiles[i], "wt", encoding="ascii") as zipfile:
                json.dump(results, zipfile)

        # record script runtime
        skm.utils.log_runtime(log[0], start_time)


rule search:
    input:
        vecfiles=rules.vectorize.output.files,
        model=join(config["model_dir"], "{fam}.pkl"),
        basis=join(config["basis_dir"], "{fam}.txt"),
        scorer=join(config["score_dir"], "{fam}.pkl")
    output:
        results=join("output", "search", "{fam}.csv")
    run:
        # load kmer basis set
        kmers = skm.io.read_output_kmers(input.basis)

        # simplify some variable names
        model_file = input.model
        score_file = input.scorer
        family = wildcards.fam

        # load model
        with open(model_file, 'rb') as mf:
            model = pickle.load(mf)

        # load scorer
        with open(score_file, 'rb') as sf:
            scorer = pickle.load(sf)

        results = list()
        for fasta in input.vecfiles:
            filename = skm.utils.split_file_ext(basename(fasta))[0]

            df = pd.read_json(fasta)
            vecs = skm.utils.to_feature_matrix(df["vector"].values)

            scores = scorer.predict(vecs, kmers)
            predictions = model.predict(scores.reshape(-1,1))
            predicted_probas = model.predict_proba(scores.reshape(-1,1))

            # display results (score, family assignment, and probability)
            df[f"{family}_score"] = scores  # scorer output
            df[family] = [True if p == 1 else False for p in predictions]
            df[f"{family}_probability"] = [p[1] for p in predicted_probas]
            df["filename"] = f"{filename}.{FILE_MAP[filename]}"

            results.append(df)

        # save full results
        results = pd.concat(results, ignore_index=True).drop(columns=['vector'])
        results.to_csv(output.results, index=False)
