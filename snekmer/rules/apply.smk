# force snakemake v6.0+ (required for modules)
from snakemake.utils import min_version

min_version("6.0")


# load snakemake modules
module process:
    snakefile:
        "process.smk"
    config:
        config


module kmerize:
    snakefile:
        "kmerize.smk"
    config:
        config


import copy
import csv

# built-in imports
import gzip
import itertools
import json
import pickle
import struct
import sys
import time
from datetime import datetime
from glob import glob
from itertools import islice, product, repeat
from multiprocessing import Pool
from os import makedirs
from os.path import basename, dirname, exists, join, split, splitext
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.csv as csv
import seaborn as sns
import sklearn
from Bio import SeqIO
from scipy import spatial
from scipy.stats import rankdata

import snekmer as skm

# Note:
# Pyarrow installed via "conda install -c conda-forge pyarrow"
# collect all fasta-like files, unzipped filenames, and basenames
input_dir = (
    "input"
    if (("input_dir" not in config) or (str(config["input_dir"]) == "None"))
    else config["input_dir"]
)
input_files = glob(join(input_dir, "*"))
# base_file = glob(join(input_dir,"base" "*"))
zipped = [fa for fa in input_files if fa.endswith(".gz")]
unzipped = [
    fa.rstrip(".gz")
    for fa, ext in product(input_files, config["input_file_exts"])
    if fa.rstrip(".gz").endswith(f".{ext}")
]
annot_files = glob(join("annotations", "*.ann"))
compare_file = glob(join("counts", "*.csv"))
confidence_file = glob(join("confidence", "*.csv"))
# map extensions to basename (basename.ext.gz -> {basename: ext})
UZ_MAP = {
    skm.utils.split_file_ext(f)[0]: skm.utils.split_file_ext(f)[1] for f in zipped
}
FA_MAP = {
    skm.utils.split_file_ext(f)[0]: skm.utils.split_file_ext(f)[1] for f in unzipped
}
# seq-annotation-scores
# get unzipped filenames
UZS = [f"{f}.{ext}" for f, ext in UZ_MAP.items()]
# isolate basenames for all files
FAS = list(FA_MAP.keys())
FAV = list(FA_MAP.values())
# parse any background files
bg_files = glob(join(input_dir, "background", "*"))
if len(bg_files) > 0:
    bg_files = [skm.utils.split_file_ext(basename(f))[0] for f in bg_files]
NON_BGS, BGS = [f for f in FAS if f not in bg_files], bg_files
# terminate with error if invalid alphabet specified
skm.alphabet.check_valid(config["alphabet"])
# define output directory (helpful for multiple runs)
out_dir = skm.io.define_output_dir(
    config["alphabet"], config["k"], nested=config["nested_output"]
)


wildcard_constraints:
    dataset=FAS,
    FAS=FAS,


# options = [(config["learnapp"]["save_apply_associations"])]
# if all((option == True or option == False) for option in options) == False:
#     sys.exit(
#         "Incorrect Value Selected. Please check if 'save_apply_associations' in in the config file under 'learnapp'. Options are 'True' or 'False'."
#     )


rule all:
    input:
        expand(join(input_dir, "{uz}"), uz=UZS),
        *[
            expand(join("output", "apply", "seq-annotation-scores-{nb}.csv"), nb=FAS)
            if config["learnapp"]["save_apply_associations"]
            else []
        ],
        expand(join("output", "apply", "kmer-summary-{nb}.csv"), nb=FAS),


use rule vectorize from kmerize with:
    input:
        fasta=lambda wildcards: join(
            input_dir, f"{wildcards.nb}.{FA_MAP[wildcards.nb]}"
        ),
    output:
        data=join("output", "vector", "{nb}.npz"),
        kmerobj=join("output", "kmerize", "{nb}.kmers"),
    log:
        join("output", "kmerize", "log", "{nb}.log"),


rule apply:
    input:
        data="output/vector/{nb}.npz",
        annotation=expand("{an}", an=annot_files),
        compare_associations=expand("{comp}", comp=compare_file),
        confidence_associations=expand("{conf}", conf=confidence_file),
    output:
        seq_ann=expand(
            join("output", "apply", "seq-annotation-scores-{nb}.csv"), nb=FAS
        )
        if config["learnapp"]["save_apply_associations"]
        else [],
        kmer_summary="output/apply/kmer-summary-{nb}.csv",
    log:
        join(out_dir, "apply", "log", "{nb}.log"),
    run:
        start_time = datetime.now()
        with open(log[0], "a") as f:
            f.write(f"start time:\t{start_time}\n")

            ##### Generate Inputs
        kmer_count_totals = pd.read_csv(
            str(input.compare_associations),
            index_col="__index_level_0__",
            header=0,
            engine="c",
        )

        kmerlist, df = skm.io.load_npz(input.data)
        kmerlist = kmerlist[0]
        seqids = df["sequence_id"]


        kmer_totals = []
        for item in kmerlist:
            kmer_totals.append(0)
        k_len = len(kmerlist[0])

        ##### Generate Kmer Counts
        seq_kmer_dict = {}
        counter = 0
        for i, seq in enumerate(seqids):
            v = df["sequence"][i]
            kmer_counts = dict()
            items = []
            for item in range(0, (len((v)) - k_len + 1)):
                items.append(v[item : (item + k_len)])
            for j in items:
                kmer_counts[j] = kmer_counts.get(j, 0) + 1
            store = []
            for i, item in enumerate(kmerlist):
                if item in kmer_counts:
                    store.append(kmer_counts[item])
                    kmer_totals[i] += kmer_counts[item]
                else:
                    store.append(0)
            seq_kmer_dict[seq] = store

            ######  Construct Kmer Counts Dataframe
        total_seqs = len(seq_kmer_dict)
        kmer_counts = pd.DataFrame(seq_kmer_dict.values())
        kmer_counts.insert(0, "Annotations", 1, True)
        kmer_totals.insert(0, total_seqs)
        kmer_counts = pd.DataFrame(
            np.insert(kmer_counts.values, 0, values=kmer_totals, axis=0)
        )
        kmer_counts.columns = ["Sequence count"] + list(kmerlist)
        kmer_counts.index = ["Totals"] + list(seq_kmer_dict.keys())

        new_associations = kmer_counts.iloc[1:, 1:].div(
            kmer_counts["Sequence count"].tolist()[1:], axis="rows"
        )

        ##### Make Kmer Counts Dataframe match Kmer Counts Totals Format
        if len(str(kmer_counts.columns.values[10])) == len(
            str(kmer_count_totals.columns.values[10])
        ):
            compare_check = True
        else:
            compare_check = False
        if compare_check == True:
            check_1 = len(new_associations.columns.values)
            check_2 = len(kmer_count_totals.columns.values)
            alphabet_initial = set(
                itertools.chain(
                    *[list(x) for x in kmer_counts.columns.values[10:check_1]]
                )
            )
            alphabet_compare = set(
                itertools.chain(
                    *[list(x) for x in kmer_count_totals.columns.values[10:check_1]]
                )
            )
            if alphabet_compare == alphabet_initial:
                compare_check = True
        if compare_check == False:
            print("Compare Check Failed. ")
            sys.exit()

        kmer_counts.drop("Totals", axis=0, inplace=True)
        kmer_counts.drop("Sequence count", axis=1, inplace=True)

        kmer_count_totals.drop("Totals", axis=0, inplace=True)
        kmer_count_totals.drop("Kmer Count", axis=1, inplace=True)
        kmer_count_totals.drop("Sequence count", axis=1, inplace=True)

        column_order = list(set(kmer_counts.columns) | set(kmer_count_totals.columns))

        kmer_counts = kmer_counts.reindex(columns=column_order, fill_value=0)
        kmer_count_totals = kmer_count_totals.reindex(
            columns=column_order, fill_value=0
        )


        #### Perform Cosine Similarity between Kmer Counts Totals and Counts and Sums DF
        cosine_df = sklearn.metrics.pairwise.cosine_similarity(
            kmer_count_totals, kmer_counts
        ).T

        kmer_count_totals = pd.DataFrame(
            cosine_df,
            columns=kmer_count_totals.index,
            index=kmer_counts.index,
        )

        ##### Write Optional Output
        if config["learnapp"]["save_apply_associations"] == True:
            out_name = (
                # "output/apply/seq-annotation-scores-" + str(input.data)[14:-4] + ".csv"
                output.seq_ann
            )
            kmer_count_totals_write = pa.Table.from_pandas(kmer_count_totals)
            csv.write_csv(kmer_count_totals_write, out_name)

            ##### Create True Output
            # Protein ID, Prediction, Score, delta, Confidence
        global_confidence_scores = pd.read_csv(str(input.confidence_associations))
        global_confidence_scores.index = global_confidence_scores[
            global_confidence_scores.columns[0]
        ]
        global_confidence_scores = global_confidence_scores.iloc[:, 1:]

        global_confidence_scores = global_confidence_scores[
            global_confidence_scores.columns[0]
        ].squeeze()

        score_rank = []
        sorted_vals = np.argsort(-kmer_count_totals.values, axis=1)[:, :2]
        for i, item in enumerate(sorted_vals):
            score_rank.append(
                (
                    kmer_count_totals[kmer_count_totals.columns[[item]]][i : i + 1]
                ).values.tolist()[0]
            )

        delta = []
        top_score = []
        for score in score_rank:
            delta.append(score[0] - score[1])
            top_score.append(score[0])

        vals = pd.DataFrame({"delta": delta})
        predictions = pd.DataFrame(kmer_count_totals.columns[sorted_vals][:, :1])
        score = pd.DataFrame(top_score)
        score.columns = ["Score"]
        predictions.columns = ["Prediction"]
        predictions = predictions.astype(str)
        vals = vals.round(decimals=2)
        vals["Confidence"] = vals["delta"].map(global_confidence_scores)

        results = pd.concat([predictions, score, vals], axis=1)
        results.index = kmer_count_totals.index

        #### Write Results
        out_name_2 = output.kmer_summary
        results.reset_index(inplace=True)
        results_write = pa.Table.from_pandas(results)
        csv.write_csv(results_write, out_name_2)

        skm.utils.log_runtime(log[0], start_time)
