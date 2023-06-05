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
from multiprocessing import pool
from os import makedirs
from os.path import basename, dirname, exists, join, split, splitext
from pathlib import path

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.csv as csv
import seaborn as sns
import sklearn
from bio import seq_io
from scipy.interpolate import interp1d
from scipy.stats import rankdata

import snekmer as skm

# note:
# pyarrow installed via "conda install -c conda-forge pyarrow"
# collect all fasta-like files, unzipped filenames, and basenames
input_dir = (
    "input"
    if (("input_dir" not in config) or (str(config["input_dir"]) == "none"))
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
uz_map = {
    skm.utils.split_file_ext(f)[0]: skm.utils.split_file_ext(f)[1] for f in zipped
}
fa_map = {
    skm.utils.split_file_ext(f)[0]: skm.utils.split_file_ext(f)[1] for f in unzipped
}

# get unzipped filenames
uzs = [f"{f}.{ext}" for f, ext in uz_map.items()]
# isolate basenames for all files
fas = list(fa_map.keys())
# parse any background files
bg_files = glob(join(input_dir, "background", "*"))
if len(bg_files) > 0:
    bg_files = [skm.utils.split_file_ext(basename(f))[0] for f in bg_files]
non_bgs, bgs = [f for f in fas if f not in bg_files], bg_files
# terminate with error if invalid alphabet specified
skm.alphabet.check_valid(config["alphabet"])
# define output directory (helpful for multiple runs)
out_dir = skm.io.define_output_dir(
    config["alphabet"], config["k"], nested=config["nested_output"]
)

options = [
    (config["learnapp"]["save_apply_associations"]),
]

if all((option == true or option == false) for option in options) == false:
    sys.exit(
        "incorrect value selected. please check a 'save_summary','save_apply_associations', or 'save_results' in the config file under 'learnapp'. options are 'true' or 'false'."
    )


# define output files to be created by snekmer
rule all:
    input:
        expand(join("output", "learn", "kmer-counts-{nb}.csv"), nb=fas),
        join("output", "learn", "kmer-counts-total.csv"),
        expand(join("output", "eval_apply", "seq-annotation-scores-{nb}.csv"), nb=fas),
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
            input_dir, f"{wildcards.nb}.{fa_map[wildcards.nb]}"
        ),
    output:
        data=join("output", "vector", "{nb}.npz"),
        kmerobj=join("output", "kmerize", "{nb}.kmers"),
    log:
        join("output", "kmerize", "log", "{nb}.log"),


# workflow to learn kmer associations
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

            ##### generate inputs
        annotation = list()
        for f in input.annotation:
            annotation.append(pd.read_table(f))
        annotations = pd.concat(annotation)
        seq_annot = {}
        seqs = annotation[0]["id"].tolist()
        anns = annotation[0]["tigrfa_ms"].tolist()
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

            ##### generate kmer counts
        k_len = len(kmerlist[0])
        seq_kmer_dict = {}
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

            # filter out non-training annotations
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

                    # construct kmer counts output
        kmer_counts = pd.data_frame(seq_kmer_dict.values())
        kmer_counts.insert(0, "annotations", annotation_counts.values(), true)
        # kmer_counts.insert(1,"kmer count",(kmer_counts[list(kmer_counts.columns[1:])].sum(axis=1).to_list()),true)
        kmer_counts_values = (
            kmer_counts[list(kmer_counts.columns[1:])].sum(axis=1).to_list()
        )
        kmer_counts.insert(1, "kmer count", kmer_counts_values, true)
        kmer_totals[0:0] = [0, total_seqs]
        colnames = ["sequence count"] + ["kmer count"] + list(kmerlist)
        kmer_counts = pd.data_frame(
            np.insert(kmer_counts.values, 0, values=(kmer_totals), axis=0)
        )
        kmer_counts.columns = colnames
        new_index = ["totals"] + list(annotation_counts.keys())
        kmer_counts.index = new_index

        #### write output
        out_name = "output/learn/kmer-counts-" + str(input.data)[14:-4] + ".csv"
        kmer_counts_out = pa.table.from_pandas(kmer_counts, preserve_index=true)
        csv.write_csv(kmer_counts_out, out_name)
        skm.utils.log_runtime(log[0], start_time)


rule merge:
    input:
        counts=expand(join("output", "learn", "kmer-counts-{nb}.csv"), nb=fas),
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
            print(kmer_counts)
            if file_num == 0:
                running_merge = kmer_counts
            elif file_num >= 1:
                running_merge = (
                    pd.concat([running_merge, kmer_counts])
                    .reset_index()
                    .groupby("__index_level_0__", sort=false)
                    .sum(min_count=1)
                ).fillna(0)

                ##### check for "base" file to merge with.
        base_check = false
        print("\n_checking for base file to merge with.\n")
        if "csv" in str(input.base_counts):
            print(
                "csv detected. matching annotations, kmers, and totals will be summed. new annotations and kmers will be added.\n"
            )
            base_check = true
        elif input.base_counts == "":
            print("no base directory detected\n")
        elif str(input.base_counts) == "base":
            print("empty base directory detected\n")
        else:
            print(
                "no file type detected. please use a .csv file in input/base directory.\n"
            )

            ##### confirm kmer counts and alphabet match base
        if base_check == true:
            base_df = pd.read_csv(
                str(input.base_counts),
                index_col="__index_level_0__",
                header=0,
                engine="pyarrow",
            )
            print("\n_base database: \n")
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
                base_check = true
            else:
                base_check = false
                print("different alphabets detected. base file not merged.")
        if base_check == true:
            print(len(str(running_merge.columns.values[1])))
            if len(str(running_merge.columns.values[1])) == len(
                str(base_df.columns.values[1])
            ):
                base_check = true
            else:
                base_check = false
                print("different kmer lengths detected. base file not merged.")

                ##### merge base file.
        if base_check == true:
            print("\n current data merged with base file.\n")
            xy = (
                pd.concat([base_df, running_merge])
                .reset_index()
                .groupby("__index_level_0__", sort=false)
                .sum(min_count=1)
            ).fillna(0)
            xy_out = pa.table.from_pandas(xy, preserve_index=true)
            csv.write_csv(xy_out, "output/learn/kmer-counts-total.csv")
            print(xy)
        else:
            print("current data merged. no detected base file to merge with. \n")
            running_merge_out = pa.table.from_pandas(running_merge, preserve_index=true)
            csv.write_csv(running_merge_out, "output/learn/kmer-counts-total.csv")

        skm.utils.log_runtime(log[0], start_time)


rule eval_apply:
    input:
        data="output/vector/{nb}.npz",
        counts="output/learn/kmer-counts-{nb}.csv",
        annotation=expand("{an}", an=annot_files),
        compare_associations=join("output", "learn", "kmer-counts-total.csv"),
    output:
        "output/eval_apply/seq-annotation-scores-{nb}.csv",
    log:
        join(out_dir, "eval_apply", "log", "{nb}.log"),
    run:
        start_time = datetime.now()
        with open(log[0], "a") as f:
            f.write(f"start time:\t{start_time}\n")

            ##### generate inputs
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
        anns = annotation[0]["tigrfa_ms"].tolist()
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

            ##### generate kmer counts
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

            ###### add known / unknown tag to mark for confidence assessment
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

        ######  construct kmer counts dataframe
        kmer_counts = pd.data_frame(seq_kmer_dict.values())
        kmer_counts.insert(0, "annotations", 1, true)
        kmer_totals.insert(0, total_seqs)
        kmer_counts = pd.data_frame(
            np.insert(kmer_counts.values, 0, values=kmer_totals, axis=0)
        )
        kmer_counts.columns = ["sequence count"] + list(kmerlist)
        kmer_counts.index = ["totals"] + list(seq_kmer_dict.keys())


        ##### make new counts data match kmer counts totals format
        if len(str(kmer_counts.columns.values[10])) == len(
            str(kmer_count_totals.columns.values[10])
        ):
            compare_check = true
        else:
            compare_check = false
        if compare_check == true:
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
                compare_check = true
            else:
                compare_check = false
        if compare_check == false:
            print("compare check failed. ")
            sys.exit()

        new_cols = set(kmer_counts.columns)
        compare_cols = set(kmer_count_totals.columns)
        add_to_compare = []
        add_to_new = []
        for val in new_cols:
            if val not in compare_cols:
                add_to_compare.append(val)
        for val in compare_cols:
            if val not in new_cols:
                add_to_new.append(val)

        kmer_count_totals = pd.concat(
            [
                kmer_count_totals,
                pd.data_frame(
                    dict.fromkeys(add_to_compare, 0), index=kmer_count_totals.index
                ),
            ],
            axis=1,
        )
        kmer_count_totals.drop(
            columns=kmer_count_totals.columns[:2], index="totals", axis=0, inplace=true
        )
        kmer_counts = pd.concat(
            [
                kmer_counts,
                pd.data_frame(dict.fromkeys(add_to_new, 0), index=kmer_counts.index),
            ],
            axis=1,
        )
        kmer_counts.drop(
            columns=kmer_counts.columns[-1:].union(kmer_counts.columns[:1]),
            index="totals",
            axis=0,
            inplace=true,
        )

        # perform cosine similarity between kmer counts totals and counts and sums df
        cosine_df = sklearn.metrics.pairwise.cosine_similarity(
            kmer_count_totals, kmer_counts
        ).t

        final_matrix_with_scores = pd.data_frame(
            cosine_df, columns=kmer_count_totals.index, index=kmer_counts.index
        )

        # write output
        out_name = (
            "output/eval_apply/seq-annotation-scores-" + str(input.data)[14:-4] + ".csv"
        )
        final_matrix_with_scores_write = pa.table.from_pandas(final_matrix_with_scores)
        csv.write_csv(final_matrix_with_scores_write, out_name)

        skm.utils.log_runtime(log[0], start_time)


rule evaluate:
    input:
        eval_apply_data=expand(
            join("output", "eval_apply", "seq-annotation-scores-{nb}.csv"), nb=fas
        ),
        base_confidence=expand("{bc}", bc=base_confidence),
    output:
        "output/eval_conf/confidence-matrix.csv",
        "output/eval_conf/global-confidence-scores.csv",
    log:
        join(out_dir, "eval_conf", "log", "conf.log"),
    run:
        start_time = datetime.now()
        with open(log[0], "a") as f:
            f.write(f"start time:\t{start_time}\n")

            #### generate input data
        print("Number of Dataframe to parse and merge: ", len(input.eval_apply_data))
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
                    tf.append("t")
                else:
                    tf.append("f")
                if "unknown" in result[i]:
                    known.append("unknown")
                else:
                    known.append("known")

            seq_ann_vals = seq_ann_scores.values
            seq_ann_vals = seq_ann_scores.values[
                np.arange(len(seq_ann_scores))[:, none],
                np.argpartition(-seq_ann_vals, np.arange(2), axis=1)[:, :2],
            ]

            diff_df = pd.data_frame(
                seq_ann_vals,
                columns=["top", "second"],
            )
            diff_df["difference"] = -(np.diff(seq_ann_vals, axis=1).round(decimals=2))
            diff_df["prediction"] = list(max_value_index)
            diff_df["actual"] = result
            diff_df["t/f"] = tf
            diff_df["known/unknown"] = known

            #### create cross_tabs - ie true/false count sums and sum within .01 intervals
            known_true_diff_df = diff_df[
                (diff_df["known/unknown"] == "known") & (diff_df["t/f"] == "t")
            ]
            known_false_diff_df = diff_df[
                (diff_df["known/unknown"] == "known") & (diff_df["t/f"] == "f")
            ]
            possible_vals = [round(x * 0.01, 2) for x in range(0, 101)]
            true_crosstab = pd.crosstab(
                known_true_diff_df.prediction, known_true_diff_df.difference
            )
            false_crosstab = pd.crosstab(
                known_false_diff_df.prediction, known_false_diff_df.difference
            )

            if j == 0:
                true_running_crosstab = true_crosstab
                false_running_crosstab = false_crosstab
            else:
                true_running_crosstab = (
                    pd.concat([true_running_crosstab, true_crosstab])
                    .reset_index()
                    .groupby("prediction", sort=false)
                    .sum(min_count=1)
                ).fillna(0)
                false_running_crosstab = (
                    pd.concat([false_running_crosstab, false_crosstab])
                    .reset_index()
                    .groupby("prediction", sort=false)
                    .sum(min_count=1)
                ).fillna(0)

            add_to_true_df = pd.data_frame(
                0,
                index=sorted(
                    set(false_running_crosstab.index) - set(true_running_crosstab.index)
                ),
                columns=true_running_crosstab.columns,
            )
            add_to_false_df = pd.data_frame(
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

            true_running_crosstab.index.names = ["prediction"]
            false_running_crosstab.index.names = ["prediction"]
            true_running_crosstab.sort_index(inplace=true)
            false_running_crosstab.sort_index(inplace=true)
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

            print("dataframes joined: ", j)
            #### generate each global cross_tab

        ratio_running_crosstab = true_running_crosstab / (
            true_running_crosstab + false_running_crosstab
        )

        true_total_dist = true_running_crosstab.sum(numeric_only=true, axis=0)
        false_total_dist = false_running_crosstab.sum(numeric_only=true, axis=0)

        print("\n_checking for base file to merge with.\n")
        print("input.base_confidence: ", input.base_confidence)
        base_conf_len = len(input.base_confidence)
        if (
            base_conf_len == 2
            and "global-true.csv" in input.base_confidence
            and "global-false.csv" in input.base_confidence
        ):
            print("two base confidence files found. will attempt to merge.")
            base_true = pd.read_csv(
                "base/confidence/global-true.csv", index_col=false, header=0
            ).set_index("difference")["0"]
            base_false = pd.read_csv(
                "base/confidence/global-false.csv", index_col=false, header=0
            ).set_index("difference")["0"]
            true_total_dist = true_total_dist + base_true
            false_total_dist = false_total_dist + base_false
        else:
            print("Base confidence files not found for merge.")

        ratio_total_dist = true_running_crosstab.sum(numeric_only=true, axis=0) / (
            true_running_crosstab.sum(numeric_only=true, axis=0)
            + false_running_crosstab.sum(numeric_only=true, axis=0)
        )

        ####interpolate for final ratio, this only will affect upper limit values if there is a decent amount of data
        ratio_total_dist = ratio_total_dist.interpolate(method="linear")

        ##### write final confidence results
        ratio_total_dist.to_csv("output/eval_conf/global-confidence-scores.csv")

        true_total_dist.to_csv("output/eval_conf/global-true.csv")

        false_total_dist.to_csv("output/eval_conf/global-false.csv")

        csv.write_csv(
            pa.table.from_pandas(ratio_running_crosstab),
            "output/eval_conf/confidence-matrix.csv",
        )

        skm.utils.log_runtime(log[0], start_time)
