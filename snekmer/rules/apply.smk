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
# Seq-Annotation-Scores
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


# check method
methods = ["score", "score_enhanced", "cosine"]
if config["learnapp"]["type"] not in methods:
    sys.exit(
        "Please select a scoring 'type' in the config file under 'learnapp'. Options are 'score', 'score_enhanced', and 'cosine'."
    )
else:
    method = config["learnapp"]["type"]
# check other learnapp config values that are T/F
options = [
    (config["learnapp"]["save_summary"]),
    (config["learnapp"]["save_apply_associations"]),
    (config["learnapp"]["save_results"]),
]
if all((option == True or option == False) for option in options) == False:
    sys.exit(
        "Incorrect Value Selected. Please check a 'save_summary','save_apply_associations', or 'save_results' in the config file under 'learnapp'. Options are 'True' or 'False'."
    )
# check summary_topN
if config["learnapp"]["save_summary"] == True:
    if type(config["learnapp"]["summary_topN"]) != int:
        sys.exit(
            "Incorrect Value Selected for summary_topN. Value must be an integer greater or equal to zero."
        )
    if config["learnapp"]["summary_topN"] < 0:
        sys.exit(
            "Incorrect Value Selected for summary_topN. Value must be an integer greater or equal to zero."
        )


rule all:
    input:
        expand(join(input_dir, "{uz}"), uz=UZS),
        # ["output/apply/Seq-Annotation-Scores-{fa}.csv".format(fa=FAS) for dataset in FAS],
        expand(join("output", "apply", "Seq-Annotation-Scores-{nb}.csv"), nb=FAS),
        expand(join("output", "apply", "kmer-summary-{nb}.csv"), nb=FAS),


use rule vectorize from kmerize with:
    input:
        fasta=lambda wildcards: join("input", f"{wildcards.nb}.{FA_MAP[wildcards.nb]}"),
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
        "output/apply/Seq-Annotation-Scores-{nb}.csv",
        "output/apply/kmer-summary-{nb}.csv",
    log:
        join(out_dir, "apply", "log", "{nb}.log"),
    run:
        start_time = datetime.now()
        with open(log[0], "a") as f:
            f.write(f"start time:\t{start_time}\n")

            ##### Generate Inputs
        Kmer_Count_Totals = pd.read_csv(
            str(input.compare_associations),
            index_col="__index_level_0__",
            header=0,
            engine="c",
        )
        df, kmerlist = skm.io.load_npz(input.data)
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
        Kmer_Counts = pd.DataFrame(seq_kmer_dict.values())
        Kmer_Counts.insert(0, "Annotations", 1, True)
        kmer_totals.insert(0, total_seqs)
        Kmer_Counts = pd.DataFrame(
            np.insert(Kmer_Counts.values, 0, values=kmer_totals, axis=0)
        )
        Kmer_Counts.columns = ["Sequence count"] + list(kmerlist)
        Kmer_Counts.index = ["Totals"] + list(seq_kmer_dict.keys())

        new_associations = Kmer_Counts.iloc[1:, 1:].div(
            Kmer_Counts["Sequence count"].tolist()[1:], axis="rows"
        )

        ##### Make Kmer Counts Dataframe match Kmer Counts Totals Format
        if len(str(Kmer_Counts.columns.values[10])) == len(
            str(Kmer_Count_Totals.columns.values[10])
        ):
            compare_check = True
        else:
            compare_check = False
        if compare_check == True:
            check_1 = len(new_associations.columns.values)
            check_2 = len(Kmer_Count_Totals.columns.values)
            alphabet_initial = set(
                itertools.chain(
                    *[list(x) for x in Kmer_Counts.columns.values[10:check_1]]
                )
            )
            alphabet_compare = set(
                itertools.chain(
                    *[list(x) for x in Kmer_Count_Totals.columns.values[10:check_1]]
                )
            )
            if alphabet_compare == alphabet_initial:
                compare_check = True
        if compare_check == False:
            print("Compare Check Failed. ")
            sys.exit()

        new_cols = set(Kmer_Counts.columns)
        compare_cols = set(Kmer_Count_Totals.columns)
        add_to_compare = []
        add_to_new = []
        for val in new_cols:
            if val not in compare_cols:
                add_to_compare.append(val)
        for val in compare_cols:
            if val not in new_cols:
                add_to_new.append(val)

        Kmer_Count_Totals = pd.concat(
            [
                Kmer_Count_Totals,
                pd.DataFrame(
                    dict.fromkeys(add_to_compare, 0), index=Kmer_Count_Totals.index
                ),
            ],
            axis=1,
        )
        Kmer_Count_Totals.drop(
            columns=Kmer_Count_Totals.columns[:2], index="Totals", axis=0, inplace=True
        )
        Kmer_Counts = pd.concat(
            [
                Kmer_Counts,
                pd.DataFrame(dict.fromkeys(add_to_new, 0), index=Kmer_Counts.index),
            ],
            axis=1,
        )
        Kmer_Counts.drop(
            columns=Kmer_Counts.columns[-1:].union(Kmer_Counts.columns[:1]),
            index="Totals",
            axis=0,
            inplace=True,
        )


        #### Perform Cosine Similarity between Kmer Counts Totals and Counts and Sums DF
        cosine_df = sklearn.metrics.pairwise.cosine_similarity(
            Kmer_Count_Totals, Kmer_Counts
        ).T
        Kmer_Count_Totals = pd.DataFrame(
            cosine_df, columns=Kmer_Count_Totals.index, index=Kmer_Counts.index
        )

        ##### Write Optional Output
        if config["learnapp"]["save_results"] == True:
            out_name = (
                "output/apply/Seq-Annotation-Scores-" + str(input.data)[14:-4] + ".csv"
            )
            Kmer_Count_Totals_write = pa.Table.from_pandas(Kmer_Count_Totals)
            csv.write_csv(Kmer_Count_Totals_write, out_name)

            ##### Create True Output
            # Protein ID, Prediction, Score, delta, Confidence
        Global_Confidence_Scores = pd.read_csv(str(input.confidence_associations))
        Global_Confidence_Scores.index = Global_Confidence_Scores[
            Global_Confidence_Scores.columns[0]
        ]
        Global_Confidence_Scores = Global_Confidence_Scores.iloc[:, 1:]
        Global_Confidence_Scores = Global_Confidence_Scores[
            Global_Confidence_Scores.columns[0]
        ].squeeze()

        Score_rank = []
        Sorted_vals = np.argsort(-Kmer_Count_Totals.values, axis=1)[:, :2]
        for i, item in enumerate(Sorted_vals):
            Score_rank.append(
                (
                    Kmer_Count_Totals[Kmer_Count_Totals.columns[[item]]][i : i + 1]
                ).values.tolist()[0]
            )

        delta = []
        Top_Score = []
        for score in Score_rank:
            delta.append(score[0] - score[1])
            Top_Score.append(score[0])

        vals = pd.DataFrame({"delta": delta})
        predictions = pd.DataFrame(Kmer_Count_Totals.columns[Sorted_vals][:, :1])
        score = pd.DataFrame(Top_Score)
        score.columns = ["Score"]
        predictions.columns = ["Prediction"]
        predictions = predictions.astype(str)
        vals = vals.round(decimals=2)
        vals["Confidence"] = vals["delta"].map(Global_Confidence_Scores)

        Results = pd.concat([predictions, score, vals], axis=1)
        Results.index = Kmer_Count_Totals.index

        #### Write Results
        out_name_2 = "output/apply/kmer-summary-" + str(input.data)[14:-4] + ".csv"
        Results.reset_index(inplace=True)
        Results_write = pa.Table.from_pandas(Results)
        csv.write_csv(Results_write, out_name_2)

        skm.utils.log_runtime(log[0], start_time)
