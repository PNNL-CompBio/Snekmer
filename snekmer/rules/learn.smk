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
from itertools import product, repeat
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
from scipy.interpolate import interp1d
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
base_counts = glob(join("base", "counts", "*.csv"))
base_confidence = glob(join("base", "confidence", "*.csv"))
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
    config["alphabet"], config["k"], nested=config["nested_output"]
)

# options = [
#     (config["learnapp"]["save_apply_associations"]),
# ]

# if all((option == True or option == False) for option in options) == False:
#     sys.exit(
#         "Incorrect Value Selected. Please check a 'save_summary','save_apply_associations', or 'save_results' in the config file under 'learnapp'. Options are 'True' or 'False'."
#     )


# define output files to be created by snekmer
rule all:
    input:
        expand(join("output", "learn", "kmer-counts-{nb}.csv"), nb=FAS),
        join("output", "learn", "kmer-counts-total.csv"),
        expand(join("output", "eval_apply", "seq-annotation-scores-{nb}.csv"), nb=FAS),
        "output/eval_conf/confidence-matrix.csv",
        "output/eval_conf/global-confidence-scores.csv",


# if any files are gzip zipped, unzip them
use rule unzip from process with:
    output:
        unzipped=join(input_dir, "{uz}"),
        zipped=join(input_dir, "zipped", "{uz}.gz"),


# build kmer count vectors for each basis set
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


# WORKFLOW to learn kmer associations
# collect all seq files and generate mega-cluster
rule learn:
    input:
        data="output/vector/{nb}.npz",
        annotation=expand("{an}", an=annot_files),
    output:
        counts="output/learn/kmer-counts-{nb}.csv",
    log:
        join(out_dir, "learn", "log", "learn-{nb}.log"),
    run:
        start_time = datetime.now()
        with open(log[0], "a") as f:
            f.write(f"start time:\t{start_time}\n")

            ##### Generate Inputs
        annotation = list()
        for f in input.annotation:
            annotation.append(pd.read_table(f))
        annotations = pd.concat(annotation)
        seq_annot = {}
        seqs = annotation[0]["id"].tolist()
        anns = annotation[0]["TIGRFAMs"].tolist()
        for i, seqid in enumerate(seqs):
            seq_annot[seqid] = anns[i]
        seqs = set(seqs)
        anns = set(anns)

        kmerlist, df = skm.io.load_npz(input.data)
        kmerlist = kmerlist[0]
        seqids = df["sequence_id"]

        kmer_totals = []
        for item in kmerlist:
            kmer_totals.append(0)

            ##### Generate Kmer Counts
        k_len = len(kmerlist[0])
        seq_kmer_dict = {}
        for i, seq in enumerate(seqids):
            v = df["sequence"][i]
            # print("V:", v)
            k_counts = dict()
            items = []
            for item in range(0, (len((v)) - k_len + 1)):
                items.append(v[item : (item + k_len)])
            for j in items:
                k_counts[j] = k_counts.get(j, 0) + 1
            store = []
            for i, item in enumerate(kmerlist):
                if item in k_counts:
                    store.append(k_counts[item])
                    kmer_totals[i] += k_counts[item]
                else:
                    store.append(0)
            seq_kmer_dict[seq] = store

            # Filter out Non-Training Annotations
            # Make sure this is happening...
        annotation_counts = {}
        total_seqs = len(seq_kmer_dict)
        for i, seqid in enumerate(list(seq_kmer_dict)):
            x = re.findall(r"\|(.*?)\|", seqid)[0]
            if x not in seqs:
                del seq_kmer_dict[seqid]
            else:
                if seq_annot[x] not in seq_kmer_dict:
                    seq_kmer_dict[seq_annot[x]] = seq_kmer_dict.pop(seqid)
                else:
                    zipped_lists = zip(
                        seq_kmer_dict.pop(seqid), seq_kmer_dict[seq_annot[x]]
                    )
                    seq_kmer_dict[seq_annot[x]] = [x + y for (x, y) in zipped_lists]
                if seq_annot[x] not in annotation_counts:
                    annotation_counts[seq_annot[x]] = 1
                else:
                    annotation_counts[seq_annot[x]] += 1

                    # Construct Kmer Counts Output
        kmer_counts = pd.DataFrame(seq_kmer_dict.values())
        kmer_counts.insert(0, "Annotations", annotation_counts.values(), True)

        kmer_counts_values = (
            kmer_counts[list(kmer_counts.columns[1:])].sum(axis=1).to_list()
        )

        kmer_counts.insert(1, "Kmer Count", kmer_counts_values, True)
        kmer_totals[0:0] = [total_seqs, sum(kmer_totals)]
        colnames = ["Sequence count"] + ["Kmer Count"] + list(kmerlist)

        kmer_counts = pd.DataFrame(
            np.insert(kmer_counts.values, 0, values=(kmer_totals), axis=0)
        )
        kmer_counts.columns = colnames
        new_index = ["Totals"] + list(annotation_counts.keys())
        kmer_counts.index = new_index

        #### Write Output
        out_name = "output/learn/kmer-counts-" + str(input.data)[14:-4] + ".csv"
        kmer_counts_out = pa.Table.from_pandas(kmer_counts, preserve_index=True)
        csv.write_csv(kmer_counts_out, out_name)
        skm.utils.log_runtime(log[0], start_time)


rule merge:
    input:
        counts=expand(join("output", "learn", "kmer-counts-{nb}.csv"), nb=FAS),
        base_counts=expand("{bf}", bf=base_counts),
    output:
        totals=join("output", "learn", "kmer-counts-total.csv"),
    log:
        join(out_dir, "learn", "log", "merge.log"),
    run:
        start_time = datetime.now()
        with open(log[0], "a") as f:
            f.write(f"start time:\t{start_time}\n")

        for file_num, f in enumerate(input.counts):
            print("database #: ", file_num, "\n")
            kmer_counts = pd.read_csv(
                str(f), index_col="__index_level_0__", header=0, engine="pyarrow"
            )
            if file_num == 0:
                running_merge = kmer_counts
            elif file_num >= 1:
                running_merge = (
                    pd.concat([running_merge, kmer_counts])
                    .reset_index()
                    .groupby("__index_level_0__", sort=False)
                    .sum(min_count=1)
                ).fillna(0)

                ##### Check for "Base" File to merge with.
        base_check = False
        print("\nChecking for base file to merge with.\n")
        if "csv" in str(input.base_counts):
            print(
                "CSV detected. Matching annotations, kmers, and totals will be summed. New annotations and kmers will be added.\n"
            )
            base_check = True
        elif input.base_counts == "":
            print("No base directory detected\n")
        elif str(input.base_counts) == "input/base":
            print("Empty base directory detected\n")
        else:
            print(
                "No file type detected. Please use a .csv file in input/base directory.\n"
            )

            ##### Confirm Kmer Counts and Alphabet Match Base
        if base_check == True:
            base_df = pd.read_csv(
                str(input.base_counts),
                index_col="__index_level_0__",
                header=0,
                engine="pyarrow",
            )
            print("\nBase Database: \n")
            print(base_df)
            check_1 = len(running_merge.columns.values)
            alphabet_initial = set(
                itertools.chain(
                    *[list(x) for x in running_merge.columns.values[3:check_1]]
                )
            )
            alphabet_base = set(
                itertools.chain(*[list(x) for x in base_df.columns.values[3:check_1]])
            )
            if alphabet_base == alphabet_initial:
                base_check = True
            else:
                base_check = False
                print("Different Alphabets Detected. Base File not merged.")
        if base_check == True:
            print(len(str(running_merge.columns.values[1])))
            if len(str(running_merge.columns.values[1])) == len(
                str(base_df.columns.values[1])
            ):
                base_check = True
            else:
                base_check = False
                print("Different kmer lengths detected. Base File not merged.")

                ##### Merge Base File.
        if base_check == True:
            print("\nMerged Database \n")
            xy = (
                pd.concat([base_df, running_merge])
                .reset_index()
                .groupby("__index_level_0__", sort=False)
                .sum(min_count=1)
            ).fillna(0)
            xy_out = pa.Table.from_pandas(xy, preserve_index=True)
            csv.write_csv(xy_out, "output/learn/kmer-counts-total.csv")
            print(xy)
        else:
            print("\Database Merged. Not merged with base file. \n")
            running_merge_out = pa.Table.from_pandas(running_merge, preserve_index=True)
            print("running_merge output:\n", running_merge)
            csv.write_csv(running_merge_out, "output/learn/kmer-counts-total.csv")

        skm.utils.log_runtime(log[0], start_time)


rule eval_apply:
    input:
        data="output/vector/{nb}.npz",
        counts="output/learn/kmer-counts-{nb}.csv",
        annotation=expand("{an}", an=annot_files),
        compare_associations=join("output", "learn", "kmer-counts-total.csv"),
    output:
        apply="output/eval_apply/seq-annotation-scores-{nb}.csv",
    log:
        join(out_dir, "eval_apply", "log", "{nb}.log"),
    run:
        start_time = datetime.now()
        with open(log[0], "a") as f:
            f.write(f"start time:\t{start_time}\n")

            ##### Generate Inputs
        annotation = list()
        kmer_count_totals = pd.read_csv(
            str(input.compare_associations),
            index_col="__index_level_0__",
            header=0,
            engine="c",
        )
        for f in input.annotation:
            annotation.append(pd.read_table(f))
        seqs = annotation[0]["id"].tolist()
        anns = annotation[0]["TIGRFAMs"].tolist()
        seq_annot = {}
        for i, seqid in enumerate(seqs):
            seq_annot[seqid] = anns[i]
        seqs = set(seqs)
        anns = set(anns)
        kmerlist, df = skm.io.load_npz(input.data)
        kmerlist = kmerlist[0]
        seqids = df["sequence_id"]
        kmer_totals = []
        for item in kmerlist:
            kmer_totals.append(0)

            ##### Generate Kmer Counts
        seq_kmer_dict = {}
        k_len = len(kmerlist[0])
        for i, seq in enumerate(seqids):
            v = df["sequence"][i]
            k_counts = dict()
            items = []
            for item in range(0, (len((v)) - k_len + 1)):
                items.append(v[item : (item + k_len)])
            for j in items:
                k_counts[j] = k_counts.get(j, 0) + 1
            store = []
            for i, item in enumerate(kmerlist):
                if item in k_counts:
                    store.append(k_counts[item])
                    kmer_totals[i] += k_counts[item]
                else:
                    store.append(0)
            seq_kmer_dict[seq] = store

            ###### ADD Known / Unknown tag to mark for confidence assessment
        annotation_counts = {}
        total_seqs = len(seq_kmer_dict)
        count = 0
        for seqid in list(seq_kmer_dict):
            x = re.findall(r"\|(.*?)\|", seqid)[0]
            if x not in seqs:
                seq_kmer_dict[(x + "_unknown_" + str(count))] = seq_kmer_dict.pop(seqid)
            else:
                seq_kmer_dict[
                    (seq_annot[x] + "_known_" + str(count))
                ] = seq_kmer_dict.pop(seqid)
            count += 1
        annotation_counts = {}
        total_seqs = len(seq_kmer_dict)

        ######  Construct Kmer Counts Dataframe
        kmer_counts = pd.DataFrame(seq_kmer_dict.values())
        kmer_counts.insert(0, "Annotations", 1, True)
        kmer_totals.insert(0, total_seqs)
        kmer_counts = pd.DataFrame(
            np.insert(kmer_counts.values, 0, values=kmer_totals, axis=0)
        )
        kmer_counts.columns = ["Sequence count"] + list(kmerlist)
        kmer_counts.index = ["Totals"] + list(seq_kmer_dict.keys())

        ##### Make New Counts Data match Kmer Counts Totals Format
        if len(str(kmer_counts.columns.values[10])) == len(
            str(kmer_count_totals.columns.values[10])
        ):
            compare_check = True
        else:
            compare_check = False
        if compare_check == True:
            check_1 = len(kmer_counts.columns.values)
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
            else:
                compare_check = False

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

        cosine_df = sklearn.metrics.pairwise.cosine_similarity(
            kmer_count_totals, kmer_counts
        ).T


        final_matrix_with_scores = pd.DataFrame(
            cosine_df, columns=kmer_count_totals.index, index=kmer_counts.index
        )

        print("New Cosine Dataframe:\n", final_matrix_with_scores)

        # Write Output
        out_name = (
            "output/eval_apply/seq-annotation-scores-" + str(input.data)[14:-4] + ".csv"
        )
        final_matrix_with_scores_write = pa.Table.from_pandas(final_matrix_with_scores)
        csv.write_csv(final_matrix_with_scores_write, out_name)

        skm.utils.log_runtime(log[0], start_time)


rule evaluate:
    input:
        eval_apply_data=expand(
            join("output", "eval_apply", "seq-annotation-scores-{nb}.csv"), nb=FAS
        ),
        base_confidence=expand("{bc}", bc=base_confidence),
    output:
        eval_cof="output/eval_conf/confidence-matrix.csv",
        eval_glob="output/eval_conf/global-confidence-scores.csv",
    log:
        join(out_dir, "eval_conf", "log", "conf.log"),
    run:
        start_time = datetime.now()
        with open(log[0], "a") as f:
            f.write(f"start time:\t{start_time}\n")

            #### Generate Input Data
        for j, f in enumerate(input.eval_apply_data):
            seq_ann_scores = pd.read_csv(
                f,
                index_col="__index_level_0__",
                header=0,
                engine="c",
            )
            max_value_index = seq_ann_scores.idxmax(axis="columns")
            result = max_value_index.keys()
            tf = list()
            known = list()

            for i, item in enumerate(list(max_value_index)):
                if item in result[i]:
                    tf.append("T")
                else:
                    tf.append("F")
                if "unknown" in result[i]:
                    known.append("Unknown")
                else:
                    known.append("Known")

            seq_ann_vals = seq_ann_scores.values
            seq_ann_vals = seq_ann_scores.values[
                np.arange(len(seq_ann_scores))[:, None],
                np.argpartition(-seq_ann_vals, np.arange(2), axis=1)[:, :2],
            ]

            diff_df = pd.DataFrame(
                seq_ann_vals,
                columns=["Top", "Second"],
            )
            diff_df["Difference"] = -(np.diff(seq_ann_vals, axis=1).round(decimals=2))
            diff_df["Prediction"] = list(max_value_index)
            diff_df["Actual"] = result
            diff_df["T/F"] = tf
            diff_df["Known/Unknown"] = known

            #### Create CrossTabs - ie True/False Count Sums and sum within .01 intervals
            known_true_diff_df = diff_df[
                (diff_df["Known/Unknown"] == "Known") & (diff_df["T/F"] == "T")
            ]
            known_false_diff_df = diff_df[
                (diff_df["Known/Unknown"] == "Known") & (diff_df["T/F"] == "F")
            ]
            possible_vals = [round(x * 0.01, 2) for x in range(0, 101)]
            true_crosstab = pd.crosstab(
                known_true_diff_df.Prediction, known_true_diff_df.Difference
            )
            false_crosstab = pd.crosstab(
                known_false_diff_df.Prediction, known_false_diff_df.Difference
            )

            if j == 0:
                true_running_crosstab = true_crosstab
                false_running_crosstab = false_crosstab
            else:
                true_running_crosstab = (
                    pd.concat([true_running_crosstab, true_crosstab])
                    .reset_index()
                    .groupby("Prediction", sort=False)
                    .sum(min_count=1)
                ).fillna(0)
                false_running_crosstab = (
                    pd.concat([false_running_crosstab, false_crosstab])
                    .reset_index()
                    .groupby("Prediction", sort=False)
                    .sum(min_count=1)
                ).fillna(0)

            add_to_true_df = pd.DataFrame(
                0,
                index=sorted(
                    set(false_running_crosstab.index) - set(true_running_crosstab.index)
                ),
                columns=true_running_crosstab.columns,
            )
            add_to_false_df = pd.DataFrame(
                0,
                index=sorted(
                    set(true_running_crosstab.index) - set(false_running_crosstab.index)
                ),
                columns=false_running_crosstab.columns,
            )

            true_running_crosstab = pd.concat([true_running_crosstab, add_to_true_df])[
                sorted(list(set(possible_vals) & set(true_running_crosstab.columns)))
            ].assign(
                **dict.fromkeys(
                    list(
                        map(
                            str,
                            sorted(
                                list(
                                    set(possible_vals)
                                    ^ set(true_running_crosstab.columns.astype(float))
                                )
                            ),
                        )
                    ),
                    0,
                )
            )
            false_running_crosstab = pd.concat(
                [false_running_crosstab, add_to_false_df]
            )[
                sorted(list(set(possible_vals) & set(false_running_crosstab.columns)))
            ].assign(
                **dict.fromkeys(
                    list(
                        map(
                            str,
                            sorted(
                                list(
                                    set(possible_vals)
                                    ^ set(false_running_crosstab.columns.astype(float))
                                )
                            ),
                        )
                    ),
                    0,
                )
            )

            true_running_crosstab.index.names = ["Prediction"]
            false_running_crosstab.index.names = ["Prediction"]
            true_running_crosstab.sort_index(inplace=True)
            false_running_crosstab.sort_index(inplace=True)
            true_running_crosstab.columns = true_running_crosstab.columns.astype(float)
            false_running_crosstab.columns = false_running_crosstab.columns.astype(
                float
            )
            true_running_crosstab = true_running_crosstab[
                sorted(true_running_crosstab.columns)
            ]
            false_running_crosstab = false_running_crosstab[
                sorted(false_running_crosstab.columns)
            ]

            print("Dataframes joined: ", j)

            #### generate each global cross_tab
        ratio_running_crosstab = true_running_crosstab / (
            true_running_crosstab + false_running_crosstab
        )

        true_total_dist = true_running_crosstab.sum(numeric_only=True, axis=0)
        false_total_dist = false_running_crosstab.sum(numeric_only=True, axis=0)
        print("true_total_dist:\n", true_total_dist)
        print("false_total_dist:\n", false_total_dist)

        print("\n Checking for base confidence files to merge with.\n")
        base_conf_len = len(input.base_confidence)
        if (
            base_conf_len == 2
            and "base/confidence/global-true.csv" in input.base_confidence
            and "base/confidence/global-false.csv" in input.base_confidence
        ):
            print("Two base confidence files found. Merging.")
            base_true = pd.read_csv(
                "base/confidence/global-true.csv", index_col=False, header=0
            ).set_index("Difference")["0"]
            base_false = pd.read_csv(
                "base/confidence/global-false.csv", index_col=False, header=0
            ).set_index("Difference")["0"]
            true_total_dist = true_total_dist + base_true
            false_total_dist = false_total_dist + base_false
        else:
            print("Base confidence files not found for merge.")

        ratio_total_dist = true_running_crosstab.sum(numeric_only=True, axis=0) / (
            true_running_crosstab.sum(numeric_only=True, axis=0)
            + false_running_crosstab.sum(numeric_only=True, axis=0)
        )

        #### interpolate for final ratio, this only will affect upper limit values if there is a decent amount of data
        ratio_total_dist = ratio_total_dist.interpolate(method="linear")

        ##### write final confidence results
        ratio_total_dist.to_csv("output/eval_conf/global-confidence-scores.csv")

        true_total_dist.to_csv("output/eval_conf/global-true.csv")

        false_total_dist.to_csv("output/eval_conf/global-false.csv")
        csv.write_csv(
            pa.Table.from_pandas(true_running_crosstab),
            "output/eval_conf/true_crosstab.csv",
        )
        csv.write_csv(
            pa.Table.from_pandas(false_running_crosstab),
            "output/eval_conf/false_crosstab.csv",
        )
        csv.write_csv(
            pa.Table.from_pandas(ratio_running_crosstab),
            "output/eval_conf/confidence-matrix.csv",
        )

        skm.utils.log_runtime(log[0], start_time)
