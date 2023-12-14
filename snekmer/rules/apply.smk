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


        class KmerCompare:
            def __init__(
                self,
                compare_associations,
                data,
                confidence_associations,
                output_seq_ann,
                output_kmer_summary,
            ):
                """
        Initialize KmerCompare with necessary file paths.

        Args:
            compare_associations (str): Path to compare associations file.
            data (str): Path to data file.
            confidence_associations (str): Path to confidence associations file.
            output_seq_ann (str): Path to output sequence annotation file.
            output_kmer_summary (str): Path to output kmer summary file.
        """
                self.compare_associations = compare_associations
                self.data = data
                self.confidence_associations = confidence_associations
                self.output_seq_ann = output_seq_ann
                self.output_kmer_summary = output_kmer_summary
                self.df = None

            def load_data(self):
                """
        Load kmer counts and sequence data from provided files.
        """
                self.kmer_count_totals = pd.read_csv(
                    str(self.compare_associations),
                    index_col="__index_level_0__",
                    header=0,
                    engine="c",
                )
                kmerlist, df = skm.io.load_npz(self.data)
                self.kmerlist = kmerlist[0]
                self.seqids = df["sequence_id"]
                self.df = df

            def generate_kmer_counts(self):
                """
        Generate k-mer counts for sequences present in the data.
        """
                self.kmer_totals = [0 for _ in self.kmerlist]
                k_len = len(self.kmerlist[0])
                self.seq_kmer_dict = {}
                for i, seq in enumerate(self.seqids):
                    v = self.df["sequence"][i]
                    kmer_counts = dict()
                    items = [
                        v[item : item + k_len] for item in range(0, len(v) - k_len + 1)
                    ]
                    for j in items:
                        kmer_counts[j] = kmer_counts.get(j, 0) + 1
                    store = [kmer_counts.get(item, 0) for item in self.kmerlist]
                    for i, item in enumerate(self.kmerlist):
                        self.kmer_totals[i] += kmer_counts.get(item, 0)
                    self.seq_kmer_dict[seq] = store

            def construct_kmer_counts_dataframe(self):
                """
        Construct a DataFrame to represent k-mer counts across sequences.
        """
                total_seqs = len(self.seq_kmer_dict)
                self.kmer_counts = pd.DataFrame(self.seq_kmer_dict.values())
                self.kmer_counts.insert(0, "Annotations", 1, True)
                self.kmer_totals.insert(0, total_seqs)
                self.kmer_counts = pd.DataFrame(
                    np.insert(
                        self.kmer_counts.values, 0, values=self.kmer_totals, axis=0
                    )
                )
                self.kmer_counts.columns = ["Sequence count"] + list(self.kmerlist)
                self.kmer_counts.index = ["Totals"] + list(self.seq_kmer_dict.keys())

            def match_kmer_counts_format(self):
                """
        Ensure that the format of the k-mer counts DataFrame matches the expected format.
        """
                if len(str(self.kmer_counts.columns.values[10])) == len(
                    str(self.kmer_count_totals.columns.values[10])
                ):
                    compare_check = True
                else:
                    compare_check = False

                if compare_check:
                    check_1 = len(self.kmer_counts.columns.values)
                    alphabet_initial = set(
                        itertools.chain(
                            *[
                                list(x)
                                for x in self.kmer_counts.columns.values[10:check_1]
                            ]
                        )
                    )
                    alphabet_compare = set(
                        itertools.chain(
                            *[
                                list(x)
                                for x in self.kmer_count_totals.columns.values[
                                    10:check_1
                                ]
                            ]
                        )
                    )
                    if alphabet_compare != alphabet_initial:
                        compare_check = False

                if not compare_check:
                    print("Compare Check Failed. ")
                    sys.exit()

                self.kmer_counts.drop("Totals", axis=0, inplace=True)
                self.kmer_counts.drop("Sequence count", axis=1, inplace=True)
                self.kmer_count_totals.drop("Totals", axis=0, inplace=True)
                self.kmer_count_totals.drop("Kmer Count", axis=1, inplace=True)
                self.kmer_count_totals.drop("Sequence count", axis=1, inplace=True)

                column_order = list(
                    set(self.kmer_counts.columns) | set(self.kmer_count_totals.columns)
                )
                self.kmer_counts = self.kmer_counts.reindex(
                    columns=column_order, fill_value=0
                )
                self.kmer_count_totals = self.kmer_count_totals.reindex(
                    columns=column_order, fill_value=0
                )

            def cosine_similarity(self):
                """
        Compute cosine similarity between kmer counts of sequences.
        """
                cosine_df = sklearn.metrics.pairwise.cosine_similarity(
                    self.kmer_count_totals, self.kmer_counts
                ).T
                self.kmer_count_totals = pd.DataFrame(
                    cosine_df,
                    columns=self.kmer_count_totals.index,
                    index=self.kmer_counts.index,
                )

            def format_and_write_output(self):
                """
        Format the results and write to specified output files.
        """
                if config["learnapp"]["save_apply_associations"]:
                    kmer_count_totals_write = pa.Table.from_pandas(
                        self.kmer_count_totals
                    )
                    csv.write_csv(kmer_count_totals_write, self.output_seq_ann)

                global_confidence_scores = pd.read_csv(
                    str(self.confidence_associations)
                )
                global_confidence_scores.index = global_confidence_scores[
                    global_confidence_scores.columns[0]
                ]
                global_confidence_scores = global_confidence_scores.iloc[:, 1:]
                global_confidence_scores = global_confidence_scores[
                    global_confidence_scores.columns[0]
                ].squeeze()

                score_rank = []
                sorted_vals = np.argsort(-self.kmer_count_totals.values, axis=1)[:, :2]
                for i, item in enumerate(sorted_vals):
                    score_rank.append(
                        (
                            self.kmer_count_totals[
                                self.kmer_count_totals.columns[[item]]
                            ][i : i + 1]
                        ).values.tolist()[0]
                    )

                delta = [score[0] - score[1] for score in score_rank]
                top_score = [score[0] for score in score_rank]

                vals = pd.DataFrame({"delta": delta})
                predictions = pd.DataFrame(
                    self.kmer_count_totals.columns[sorted_vals][:, :1]
                )
                score = pd.DataFrame(top_score)
                score.columns = ["Score"]
                predictions.columns = ["Prediction"]
                predictions = predictions.astype(str)
                vals = vals.round(decimals=2)
                vals["Confidence"] = vals["delta"].map(global_confidence_scores)

                results = pd.concat([predictions, score, vals], axis=1)
                results.index = self.kmer_count_totals.index

                results.reset_index(inplace=True)
                results_write = pa.Table.from_pandas(results)
                csv.write_csv(results_write, self.output_kmer_summary)

            def execute_all(self):
                """
        Execute the entire sequence of operations in the KmerCompare process.
        """
                self.load_data()
                self.generate_kmer_counts()
                self.construct_kmer_counts_dataframe()
                self.match_kmer_counts_format()
                self.cosine_similarity()
                self.format_and_write_output()


        apply = KmerCompare(
            input.compare_associations,
            input.data,
            input.confidence_associations,
            output.seq_ann,
            output.kmer_summary,
        )
        apply.execute_all()

        skm.utils.log_runtime(log[0], start_time)
