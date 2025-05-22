# ---------------------------------------------------------
# Imports
# ---------------------------------------------------------


print("Starting imports")

import pickle
from datetime import datetime
from os import makedirs
from os.path import exists, join

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import snekmer as skm
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split, StratifiedKFold
from sklearn.preprocessing import LabelEncoder

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.csv as csv
import seaborn as sns
import sklearn
from Bio import SeqIO
from scipy.interpolate import interp1d
from scipy.stats import rankdata
import random
import snekmer as skm
from scipy.ndimage import gaussian_filter1d
import sklearn.metrics.pairwise
import itertools
import os
import shutil
import sys
import re
# ---------------------------------------------------------
# Files and Parameters
# ---------------------------------------------------------

print("Starting files params")

config = snakemake.config

# change matplotlib backend to non-interactive
plt.switch_backend("Agg")

outDir = skm.io.define_output_dir(
    config["alphabet"], config["k"], nested=config["nested_output"]
)

# ---------------------------------------------------------
# Run script
# ---------------------------------------------------------


class KmerCompare:
    def __init__(
        self,
        compare_associations,
        data,
        confidence_associations,
        decoy_stats,
        annotation,
        output_seq_ann,
        output_kmer_summary,
        selection_type,
        threshold_type,
    ):
        """
Initialize KmerCompare with necessary file paths.

Args:
    compare_associations (str): Path to compare associations file.
    data (str): Path to data file.
    confidence_associations (str): Path to confidence associations file.
    decoy_stats (str): Path to decoy stats file.
    annotation (str): Path to annotation file.
    output_seq_ann (str): Path to output sequence annotation file.
    output_kmer_summary (str): Path to output kmer summary file.
    selection_type (str): Method selection based on config.
    threshold_type (str): Threshold type from config.
"""
        self.compare_associations = compare_associations
        self.data = data
        self.confidence_associations = confidence_associations
        self.annotation = annotation
        self.decoy_stats = decoy_stats
        self.output_seq_ann = output_seq_ann
        self.output_kmer_summary = output_kmer_summary
        self.selection_type = selection_type
        self.threshold_type = threshold_type
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

        # Method 0: Hit Hit No Threshold

    def select_top_no_threshold(self):
        self.selected_values = {}
        for row_id, row in self.kmer_count_totals.iterrows():
            if not row.empty:
                sorted_row = row.sort_values(ascending=False)
                top_value = sorted_row.iloc[0]
                top_family = sorted_row.index[0]
                if len(sorted_row) > 1:
                    second_value = sorted_row.iloc[1]
                    delta = top_value - second_value
                else:
                    delta = top_value
                self.selected_values[row_id] = (top_family, top_value, delta)
            else:
                self.selected_values[row_id] = (None, None, None)

                # Method 1: Top Hit Above Threshold

    def select_top_above_threshold(self):
        self.decoy_df = pd.read_csv(
            str(self.decoy_stats),
            header=0,
            engine="c",
        )
        threshold_dict = dict(
            zip(self.decoy_df.Family, self.decoy_df[self.threshold_type])
        )

        self.selected_values = {}
        filtered_out_count = 0

        for row_id, row in self.kmer_count_totals.iterrows():
            threshold_values = row.index.map(threshold_dict.get)
            threshold_series = pd.Series(threshold_values, index=row.index)

            row_values = row[row > threshold_series]

            if not row_values.empty:
                sorted_row = row_values.sort_values(ascending=False)
                top_value = sorted_row.iloc[0]
                top_family = sorted_row.index[0]
                if len(sorted_row) > 1:
                    second_value = sorted_row.iloc[1]
                    delta = top_value - second_value
                else:
                    delta = top_value - threshold_dict.get(top_family, 0)
                self.selected_values[row_id] = (top_family, top_value, delta)
            else:
                self.selected_values[row_id] = (None, None, None)
                filtered_out_count += 1

                # Method 2: Greatest Distance
                # Note, Delta is calculated different from Top Two Hit method - and as such, is not comparable in the same way.
                # Note, Score is also maybe calculated differently.  Actually I think it might still be cosine score.  CONFIRM THIS.

    def select_by_greatest_distance(self):
        self.decoy_df = pd.read_csv(
            str(self.decoy_stats),
            header=0,
            engine="c",
        )
        threshold_dict = dict(
            zip(self.decoy_df.Family, self.decoy_df[self.threshold_type])
        )

        self.selected_values = {}
        filtered_out_count = 0

        for row_id, row in self.kmer_count_totals.iterrows():
            distances = row - row.index.map(threshold_dict.get)
            positive_distances = distances[distances > 0]

            if not positive_distances.empty:
                greatest_distance_family = positive_distances.idxmax()
                greatest_distance_value = row[greatest_distance_family]
                delta = positive_distances.max()
                self.selected_values[row_id] = (
                    greatest_distance_family,
                    greatest_distance_value,
                    delta,
                )
            else:
                self.selected_values[row_id] = (None, None, None)
                filtered_out_count += 1

                # Method 3: Balanced Distance
                # Note, Delta is calculated different from Top Two Hit method - and as such, is not comparable in the same way.
                # Note, Score is also now calculated differently. It is not the Cosine Similarity Score, it is weighted

    def select_by_balanced_distance(self, weight_top=0.5, weight_distance=0.5):
        self.decoy_df = pd.read_csv(
            str(self.decoy_stats),
            header=0,
            engine="c",
        )
        threshold_dict = dict(
            zip(self.decoy_df.Family, self.decoy_df[self.threshold_type])
        )

        self.selected_values = {}
        filtered_out_count = 0

        for row_id, row in self.kmer_count_totals.iterrows():
            threshold_values = row.index.map(threshold_dict.get)
            threshold_series = pd.Series(threshold_values, index=row.index)

            row_values_above_threshold = row[row > threshold_series]
            distances = row - threshold_series
            positive_distances = distances[distances > 0]

            candidates = {}

            if not row_values_above_threshold.empty:
                top_value = row_values_above_threshold.max()
                top_family = row_values_above_threshold.idxmax()
                top_threshold = threshold_dict.get(top_family, 0)
                candidates[top_family] = {
                    "value": top_value,
                    "delta": top_value - top_threshold,
                    "score": (top_value * weight_top)
                    + ((top_value - top_threshold) * weight_distance),
                }

            if not positive_distances.empty:
                greatest_distance_family = positive_distances.idxmax()
                greatest_distance_value = row[greatest_distance_family]
                greatest_distance_threshold = threshold_dict.get(
                    greatest_distance_family, 0
                )
                candidates[greatest_distance_family] = {
                    "value": greatest_distance_value,
                    "delta": positive_distances.max(),
                    "score": (greatest_distance_value * weight_top)
                    + (positive_distances.max() * weight_distance),
                }

            if candidates:
                best_candidate = max(
                    candidates.items(), key=lambda x: x[1]["score"]
                )
                self.selected_values[row_id] = (
                    best_candidate[0],
                    best_candidate[1]["value"],
                    best_candidate[1]["delta"],
                )
            else:
                self.selected_values[row_id] = (None, None, None)
                filtered_out_count += 1

    def format_and_write_output(self):
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

        results_list = []
        for row_id in self.kmer_count_totals.index:
            if row_id in self.selected_values:
                prediction, score, delta = self.selected_values[row_id]
                if delta is None:
                    delta = 0
            else:
                prediction, score, delta = None, None, 0
            results_list.append(
                {
                    "Sequence": row_id,
                    "Prediction": prediction,
                    "Score": score,
                    "delta": round(delta, 2),
                }
            )

        results = pd.DataFrame(results_list)
        results.set_index("Sequence", inplace=True)

        results["Confidence"] = results["delta"].map(global_confidence_scores)

        results.reset_index(inplace=True)
        results_write = pa.Table.from_pandas(results)
        csv.write_csv(results_write, self.output_kmer_summary)

    def execute_all(self):
        self.load_data()
        self.generate_kmer_counts()
        self.construct_kmer_counts_dataframe()
        self.match_kmer_counts_format()
        self.cosine_similarity()
        # Select method based on selection_type
        if self.selection_type == "top_hit":
            if self.threshold_type is None:
                self.select_top_no_threshold()
            else:
                self.select_top_above_threshold()
        elif self.selection_type == "greatest_distance":
            self.select_by_greatest_distance()
        elif self.selection_type == "combined_distance":
            weight_top = config["learnapp"].get("weight_top", 0.5)
            weight_distance = config["learnapp"].get("weight_distance", 0.5)
            self.select_by_balanced_distance(weight_top, weight_distance)
        else:
            raise ValueError(f"Invalid selection_type: {self.selection_type}")
        self.format_and_write_output()


apply = KmerCompare(
    snakemake.input.compare_associations,
    snakemake.input.data,
    snakemake.input.confidence_associations,
    snakemake.input.decoy_stats,
    snakemake.input.annotation,
    snakemake.output.seq_ann,
    snakemake.output.kmer_summary,
    snakemake.params.selection_type,
    snakemake.params.threshold_type,
)
apply.execute_all()
