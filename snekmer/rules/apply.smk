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

import matplotlib.pyplot as plt
from matplotlib_venn import venn3

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
# confidence_file = glob(join("confidence", "*.csv"))
decoy_stats_file = glob(join("stats", "*.csv"))
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

threshold_type = (
    config["learnapp"]["threshold"]
)
selection_type = (
    config["learnapp"]["selection"]
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
        # confidence_associations=expand("{conf}", conf=confidence_file), # modifed
        decoy_stats = expand("{decoy}", decoy=decoy_stats_file),
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
                decoy_stats,  # modified 
                annotation,
                output_seq_ann,
                output_kmer_summary,
            ):
                """
        Initialize KmerCompare with necessary file paths.

        Args:
            compare_associations (str): Path to compare associations file.
            data (str): Path to data file.
            decoy_stats (str): Path to decoy_stats associations file. # modifed
            output_seq_ann (str): Path to output sequence annotation file.
            output_kmer_summary (str): Path to output kmer summary file.
        """
                self.compare_associations = compare_associations
                self.data = data
                # self.confidence_associations = confidence_associations   # modifed 
                self.annotation = annotation
                self.decoy_stats = decoy_stats
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

        #     def format_and_write_output(self):
        #         """
        # Format the results and write to specified output files.
        # """
        #         if config["learnapp"]["save_apply_associations"]:
        #             kmer_count_totals_write = pa.Table.from_pandas(
        #                 self.kmer_count_totals
        #             )
        #             csv.write_csv(kmer_count_totals_write, self.output_seq_ann)

        #         global_confidence_scores = pd.read_csv(
        #             str(self.confidence_associations)  # old method
        #         )
        #         global_confidence_scores.index = global_confidence_scores[
        #             global_confidence_scores.columns[0]
        #         ]
        #         global_confidence_scores = global_confidence_scores.iloc[:, 1:]
        #         global_confidence_scores = global_confidence_scores[
        #             global_confidence_scores.columns[0]
        #         ].squeeze()

        #         score_rank = []
        #         sorted_vals = np.argsort(-self.kmer_count_totals.values, axis=1)[:, :2]
        #         for i, item in enumerate(sorted_vals):
        #             score_rank.append(
        #                 (
        #                     self.kmer_count_totals[
        #                         self.kmer_count_totals.columns[[item]]
        #                     ][i : i + 1]
        #                 ).values.tolist()[0]
        #             )

        #         delta = [score[0] - score[1] for score in score_rank]
        #         top_score = [score[0] for score in score_rank]

        #         vals = pd.DataFrame({"delta": delta})
        #         predictions = pd.DataFrame(
        #             self.kmer_count_totals.columns[sorted_vals][:, :1]
        #         )
        #         score = pd.DataFrame(top_score)
        #         score.columns = ["Score"]
        #         predictions.columns = ["Prediction"]
        #         predictions = predictions.astype(str)
        #         vals = vals.round(decimals=2)
        #         vals["Confidence"] = vals["delta"].map(global_confidence_scores)

        #         results = pd.concat([predictions, score, vals], axis=1)
        #         results.index = self.kmer_count_totals.index

        #         results.reset_index(inplace=True)
        #         results_write = pa.Table.from_pandas(results)
        #         csv.write_csv(results_write, self.output_kmer_summary)


            #method_0
            def select_top_no_threshold(self):
                # Initialize an empty dictionary to store the top values
                self.top_values_no_thresh = {}
                
                # Iterate over each row in the kmer_count_totals DataFrame
                for row_id, row in self.kmer_count_totals.iterrows():
                    # Filter values (no threshold applied)
                    row_values = row
                    
                    # Get the maximum value and its corresponding family
                    if not row_values.empty:
                        top_value = row_values.max()  # Select the highest value
                        top_family = row_values.idxmax()  # Get the corresponding family (column name)
                        self.top_values_no_thresh[row_id] = (top_family, top_value)
                
                # Print a sample of top values for debugging
                print(f"Top values: {dict(itertools.islice(self.top_values_no_thresh.items(), 10))}")
                print(f"len(self.top_values_no_thresh): {len(self.top_values_no_thresh)}")


            #Method 1
            def select_top_above_threshold(self):
                self.decoy_df = pd.read_csv(
                    str(self.decoy_stats),
                    header=0,
                    engine="c",
                )
                threshold_dict = dict(zip(self.decoy_df.Family, self.decoy_df[threshold_type]))
                
                self.top_values = {}
                filtered_out_count = 0
                top_filtered_count = 0

                # for row_id, row in self.kmer_count_totals.iterrows():
                #     # Filter values above threshold
                #     row_values = row[row > row.index.map(threshold_dict.get)]
                #     if not row_values.empty:
                #         top_value = row_values.max()  # Select the highest value above the threshold
                #         top_family = row_values.idxmax()  # Get the corresponding family (column name)
                #         self.top_values[row_id] = (top_family, top_value)
                #     else:
                #         filtered_out_count += 1

                for row_id, row in self.kmer_count_totals.iterrows():
                    # Get threshold values for the current row family
                    threshold_values = row.index.map(threshold_dict.get)
                    
                    # Retain values above the threshold; values below the threshold and NaNs are excluded
                    row_values = row[row > threshold_values]

                    if not row_values.empty:
                        top_value = row_values.max()  # Select the highest value above the threshold
                        top_family = row_values.idxmax()  # Get the corresponding family (column name)

                        # Check if the top_value is the highest in the entire row
                        if top_value != row.max():  # Compare with the maximum value of the original row
                            # print(f"Family: {top_family}, Threshold: {threshold_values[row.index.get_loc(top_family)]}, Value: {top_value}")
                            top_filtered_count +=1
                        self.top_values[row_id] = (top_family, top_value)
                    else:
                        filtered_out_count += 1
                
                print(f"top values: {dict(itertools.islice(self.top_values.items(), 10))}")
                print(f"len(self.top_values): {len(self.top_values)}")
                print(f"Full Rows filtered by threshold: {filtered_out_count}")
                print(f"Top Values filtered by threshold: {top_filtered_count}")




            # Method 2: Select value with the greatest distance from its threshold
            def select_by_greatest_distance(self):
                self.decoy_df = pd.read_csv(
                    str(self.decoy_stats),
                    header=0,
                    engine="c",
                )
                threshold_dict = dict(zip(self.decoy_df.Family, self.decoy_df[threshold_type]))
                
                self.distance_values = {}
                filtered_out_count = 0
                
                for row_id, row in self.kmer_count_totals.iterrows():
                    # Calculate distance from threshold
                    distances = row - row.index.map(threshold_dict.get)
                    positive_distances = distances[distances > 0]  # Only consider values above their thresholds
                    
                    if not positive_distances.empty:
                        greatest_distance_family = positive_distances.idxmax()  # Get the family with the greatest distance
                        greatest_distance_value = row[greatest_distance_family]  # Get the corresponding value
                        self.distance_values[row_id] = (greatest_distance_family, greatest_distance_value)
                    else:
                        filtered_out_count += 1

                print(f"distance values: {dict(itertools.islice(self.distance_values.items(), 10))}")
                print(f"Distance vals filtered by threshold: {filtered_out_count}")


            # #Method 3: Relative Distance Above Threshold
            # def select_by_greatest_relative_distance(self):
            #     # Load decoy statistics
            #     self.decoy_df = pd.read_csv(
            #         str(self.decoy_stats),
            #         header=0,
            #         engine="c",
            #     )
            #     threshold_dict = dict(zip(self.decoy_df.Family, self.decoy_df[threshold_type]))
                
            #     self.relative_distance_values = {}
            #     filtered_out_count = 0
                
            #     for row_id, row in self.kmer_count_totals.iterrows():
            #         # Calculate relative distance from threshold
            #         relative_distances = (row - row.index.map(threshold_dict.get)) / row.index.map(threshold_dict.get)
            #         positive_relative_distances = relative_distances[relative_distances > 0]  # Only consider positive relative distances (i.e., values above their thresholds)
                    
            #         if not positive_relative_distances.empty:
            #             greatest_relative_distance_family = positive_relative_distances.idxmax()  # Get the family with the greatest relative distance
            #             greatest_relative_distance_value = row[greatest_relative_distance_family]  # Get the corresponding value
            #             self.relative_distance_values[row_id] = (greatest_relative_distance_family, greatest_relative_distance_value)
            #         else:
            #             filtered_out_count += 1

            #     # Print out top results for inspection
            #     print(f"relative distance values: {dict(itertools.islice(self.relative_distance_values.items(), 10))}")
            #     print(f"Relative distance vals filtered by threshold: {filtered_out_count}")

            # #Method 4: 
            # def select_by_combined_distance(self, weight_absolute=0.5, weight_relative=0.5):
            #     """
            #     Selects values based on a balance between absolute and relative distance from thresholds.
            #     Parameters:
            #         - weight_absolute: weight given to absolute distance
            #         - weight_relative: weight given to relative distance
            #     """
            #     # Load decoy statistics
            #     self.decoy_df = pd.read_csv(
            #         str(self.decoy_stats),
            #         header=0,
            #         engine="c",
            #     )
            #     threshold_dict = dict(zip(self.decoy_df.Family, self.decoy_df[threshold_type]))

            #     self.combined_distance_values = {}
            #     filtered_out_count = 0

            #     for row_id, row in self.kmer_count_totals.iterrows():
            #         # Calculate absolute distance
            #         absolute_distances = row - row.index.map(threshold_dict.get)
            #         positive_absolute_distances = absolute_distances[absolute_distances > 0]  # Only consider positive distances
                    
            #         # Calculate relative distance
            #         relative_distances = (row - row.index.map(threshold_dict.get)) / row.index.map(threshold_dict.get)
            #         positive_relative_distances = relative_distances[relative_distances > 0]  # Only consider positive relative distances

            #         # Combine absolute and relative distances with given weights
            #         combined_distances = (positive_absolute_distances * weight_absolute) + (positive_relative_distances * weight_relative)

            #         if not combined_distances.empty:
            #             greatest_combined_distance_family = combined_distances.idxmax()  # Get the family with the greatest combined distance
            #             greatest_combined_distance_value = row[greatest_combined_distance_family]  # Get the corresponding value
            #             self.combined_distance_values[row_id] = (greatest_combined_distance_family, greatest_combined_distance_value)
            #         else:
            #             filtered_out_count += 1

            #     # Print top results for inspection
            #     print(f"Combined distance values: {dict(itertools.islice(self.combined_distance_values.items(), 10))}")
            #     print(f"Values filtered by threshold: {filtered_out_count}")

            #Method 5
            def select_by_combined_top_and_distance(self, weight_top=0.5, weight_distance=0.5):
                """
                Selects values based on a balance between the top value above the threshold and the value with the greatest distance from the threshold.
                Parameters:
                    - weight_top: weight given to selecting the top value above the threshold
                    - weight_distance: weight given to selecting the value with the greatest distance from the threshold
                """
                # Load decoy statistics
                self.decoy_df = pd.read_csv(
                    str(self.decoy_stats),
                    header=0,
                    engine="c",
                )
                threshold_dict = dict(zip(self.decoy_df.Family, self.decoy_df[threshold_type]))

                self.combined_top_distance_values = {}
                filtered_out_count = 0

                for row_id, row in self.kmer_count_totals.iterrows():
                    # Get the threshold values for the current row
                    threshold_values = row.index.map(threshold_dict.get)
                    # Convert threshold_values to a pandas Series with the same index as row
                    threshold_values = pd.Series(threshold_values, index=row.index)

                    # Step 1: Get top value above threshold
                    row_values_above_threshold = row[row > threshold_values]
                    if not row_values_above_threshold.empty:
                        top_value = row_values_above_threshold.max()  # Select the highest value above the threshold
                        top_family = row_values_above_threshold.idxmax()  # Get the family with the top value above threshold
                    else:
                        top_value, top_family = None, None

                    # Step 2: Calculate absolute distance from threshold
                    distances = row - threshold_values
                    positive_distances = distances[distances > 0]  # Only consider positive distances
                    if not positive_distances.empty:
                        greatest_distance_family = positive_distances.idxmax()  # Get the family with the greatest distance
                        greatest_distance_value = row[greatest_distance_family]  # Get the corresponding value
                    else:
                        greatest_distance_family, greatest_distance_value = None, None

                    # Step 3: Combine both strategies using weighted scores
                    if top_value is not None and greatest_distance_value is not None:
                        # Calculate a combined score for each family based on both the cosine similarity (value) and distance from threshold
                        top_threshold = threshold_dict.get(top_family)
                        greatest_distance_threshold = threshold_dict.get(greatest_distance_family)

                        # Ensure thresholds are not None
                        if top_threshold is None or greatest_distance_threshold is None:
                            filtered_out_count += 1
                            continue

                        top_score = (top_value * weight_top) + ((top_value - top_threshold) * weight_distance)
                        greatest_distance_score = (greatest_distance_value * weight_top) + ((greatest_distance_value - greatest_distance_threshold) * weight_distance)

                        # Select the family with the highest combined score
                        if top_score >= greatest_distance_score:
                            combined_family = top_family
                            combined_value = top_value
                        else:
                            combined_family = greatest_distance_family
                            combined_value = greatest_distance_value

                    elif top_value is not None:
                        # If only the top value is available, use the top family
                        combined_family = top_family
                        combined_value = top_value
                    elif greatest_distance_value is not None:
                        # If only the greatest distance value is available, use the greatest distance family
                        combined_family = greatest_distance_family
                        combined_value = greatest_distance_value
                    else:
                        filtered_out_count += 1
                        continue

                    # Store the selected family and value as a tuple
                    self.combined_top_distance_values[row_id] = (combined_family, combined_value)

                # Print top results for inspection
                print(f"Combined top and distance values: {dict(itertools.islice(self.combined_top_distance_values.items(), 10))}")
                print(f"Values filtered by threshold: {filtered_out_count}")


            # Self Validation Script
            def calculate_accuracy_greatest_distance(self):

                def _extract_id(full_id):
                    """Extracts the ID from a string, splitting on '|' and returning the second part if it exists."""
                    return full_id.split('|')[1] if '|' in full_id else full_id

                def _evaluate_method(annotation_dict, prediction_dict, method_name):
                    total_known = 0
                    total_unknown = 0
                    total = 0
                    true = 0
                    false = 0
                    
                    for k, v in prediction_dict.items():
                        fam_id = _extract_id(k)
                        predicted_family = v[0] if isinstance(v, tuple) else v  # Adjusted to handle both formats
                        if fam_id in annotation_dict:
                            if annotation_dict[fam_id] == predicted_family:
                                true += 1
                            else:
                                false += 1
                            total_known += 1
                        else:
                            total_unknown += 1
                        total += 1
                    
                    # Calculate accuracy
                    accuracy = true / (true + false) if (true + false) > 0 else 0

                    print(f"---- {method_name} Summary ----")
                    print(f"True: {true}")
                    print(f"False: {false}")
                    print(f"Unknown: {total_unknown}")
                    print(f"Total known: {total_known}")
                    print(f"Total Items: {len(prediction_dict)}")
                    print(f"Accuracy: {accuracy:.4f}")

                    return {
                        'method': method_name,
                        'true': true,
                        'false': false,
                        'unknown': total_unknown,
                        'known': total_known,
                        'accuracy': accuracy
                    }

                # Load the annotations
                annotation_df = pd.read_csv(
                    str(self.annotation),
                    header=0,
                    engine="c",
                    sep="\t"
                )

                annotation_dict = dict(zip(annotation_df['id'], annotation_df['TIGRFAMs']))

                # Evaluate existing methods
                distance_eval = _evaluate_method(annotation_dict, self.distance_values, "Greatest Distance")
                top_eval = _evaluate_method(annotation_dict, self.top_values, "Top Values")
                no_thresh_eval = _evaluate_method(annotation_dict, self.top_values_no_thresh, "Top Values No Threshold")
                rel_dist_eval = _evaluate_method(annotation_dict, self.relative_distance_values, "Relative Distance")

                # Create a summary table
                summary_table = pd.DataFrame([distance_eval, top_eval, no_thresh_eval, rel_dist_eval])

                weight_combinations = [
                    (0.1, 0.9),
                    (0.2, 0.8),
                    (0.3, 0.7),
                    (0.4, 0.6),
                    (0.5, 0.5),
                    (0.6, 0.4),
                    (0.7, 0.3),
                    (0.8, 0.2),
                    (0.9, 0.1),
                    (0.99, 0.01),
                    (1.0, 0.0)  # Pure absolute distance for comparison
                ]

                # # Loop over weight combinations for combined distance method
                # for i, (weight_absolute, weight_relative) in enumerate(weight_combinations):
                #     # Use the combined distance method with the current weight combination
                #     self.select_by_combined_distance(weight_absolute=weight_absolute, weight_relative=weight_relative)

                #     # Evaluate the combined distance method
                #     combined_eval = _evaluate_method(
                #         annotation_dict, 
                #         self.combined_distance_values, 
                #         f"Combined Distance ({weight_absolute},{weight_relative})"
                #     )
                    
                #     # Add the evaluation to the summary table
                #     summary_table = summary_table.append(combined_eval, ignore_index=True)

                # Loop over weight combinations for combined top and distance method
                for i, (weight_top, weight_distance) in enumerate(weight_combinations):
                    # Use the combined top and distance method with the current weight combination
                    self.select_by_combined_top_and_distance(weight_top=weight_top, weight_distance=weight_distance)

                    # Evaluate the combined top and distance method
                    combined_top_distance_eval = _evaluate_method(
                        annotation_dict, 
                        self.combined_top_distance_values, 
                        f"Combined Top and Distance ({weight_top},{weight_distance})"
                    )
                    
                    # Add the evaluation to the summary table
                    summary_table = summary_table.append(combined_top_distance_eval, ignore_index=True)

                print("\n--- Summary Table ---")
                print(summary_table)
                                                            


                # # Prepare sets for true and false predictions
                # true_distance_set = set()
                # false_distance_set = set()

                # true_top_set = set()
                # false_top_set = set()

                # true_rel_dist_set = set()
                # false_rel_dist_set = set()

                # true_no_thresh_set = set()
                # false_no_thresh_set = set()

                # # Categorize predictions as true or false for each method
                # for k, v in self.distance_values.items():
                #     fam_id = _extract_id(k)
                #     if fam_id in annotation_dict:
                #         if annotation_dict[fam_id] == v[0]:
                #             true_distance_set.add(k)
                #         else:
                #             false_distance_set.add(k)

                # for k, v in self.top_values.items():
                #     fam_id = _extract_id(k)
                #     if fam_id in annotation_dict:
                #         if annotation_dict[fam_id] == v[0]:
                #             true_top_set.add(k)
                #         else:
                #             false_top_set.add(k)

                # for k, v in self.relative_distance_values.items():
                #     fam_id = _extract_id(k)
                #     if fam_id in annotation_dict:
                #         if annotation_dict[fam_id] == v[0]:
                #             true_rel_dist_set.add(k)
                #         else:
                #             false_rel_dist_set.add(k)

                # # Categorize predictions for top_values_no_thresh
                # for k, v in self.top_values_no_thresh.items():
                #     fam_id = _extract_id(k)
                #     if fam_id in annotation_dict:
                #         if annotation_dict[fam_id] == v[0]:
                #             true_no_thresh_set.add(k)
                #         else:
                #             false_no_thresh_set.add(k)

                # # Use the Agg backend for headless operation
                # plt.switch_backend('Agg')

                # Plot Venn diagrams for TRUE predictions overlap
                # Helper function to print colors used in Venn diagram
                # Helper function to print colors of main labels in the Venn diagram

                # def print_label_colors(venn, labels):
                #     """Print the colors of the main labels in the Venn diagram."""
                #     for subset, label in zip(['100', '010', '001'], labels):
                #         if venn.get_patch_by_id(subset) is not None:
                #             color = venn.get_patch_by_id(subset).get_facecolor()[:3]  # Extract RGB color
                #             print(f"{label} Color: {color}")

                # # For TRUE predictions overlap
                # plt.figure(figsize=(8, 8))
                # venn_true_orig = venn3([true_distance_set, true_top_set, true_rel_dist_set], 
                #                         ('Absolute Distance', 'Top Values', 'Relative Distance'))
                # plt.title('Overlap of Correct Predictions: Absolute Distance vs Top Values vs Relative Distance')
                # print("Colors for 'Absolute Distance', 'Top Values', 'Relative Distance'")
                # print_label_colors(venn_true_orig, ['Absolute Distance', 'Top Values', 'Relative Distance'])
                # plt.savefig('venn_diagram_true_original_three.png')
                # plt.close()

                # # Replace Relative Distance with Top Values No Thresh
                # plt.figure(figsize=(8, 8))
                # venn_true_no_thresh = venn3([true_distance_set, true_top_set, true_no_thresh_set], 
                #                             ('Absolute Distance', 'Top Values', 'Top Values No Thresh'))
                # plt.title('Overlap of Correct Predictions: Absolute Distance vs Top Values vs Top Values (No Threshold)')
                # print("Colors for 'Absolute Distance', 'Top Values', 'Top Values (No Thresh)'")
                # print_label_colors(venn_true_no_thresh, ['Absolute Distance', 'Top Values', 'Top Values (No Thresh)'])
                # plt.savefig('venn_diagram_true_no_thresh.png')
                # plt.close()

                # # For FALSE predictions overlap
                # plt.figure(figsize=(8, 8))
                # venn_false_orig = venn3([false_distance_set, false_top_set, false_rel_dist_set], 
                #                         ('Absolute Distance', 'Top Values', 'Relative Distance'))
                # plt.title('Overlap of Incorrect Predictions: Absolute Distance vs Top Values vs Relative Distance')
                # print("Colors for Incorrect Predictions: 'Absolute Distance', 'Top Values', 'Relative Distance'")
                # print_label_colors(venn_false_orig, ['Absolute Distance', 'Top Values', 'Relative Distance'])
                # plt.savefig('venn_diagram_false_original_three.png')
                # plt.close()

                # Replace Relative Distance with Top Values No Thresh
                # plt.figure(figsize=(8, 8))
                # venn_false_no_thresh = venn3([false_distance_set, false_top_set, false_no_thresh_set], 
                #                             ('Absolute Distance', 'Top Values', 'Top Values No Thresh'))
                # plt.title('Overlap of Incorrect Predictions: Absolute Distance vs Top Values vs Top Values (No Threshold)')
                # print("Colors for Incorrect Predictions: 'Absolute Distance', 'Top Values', 'Top Values (No Thresh)'")
                # print_label_colors(venn_false_no_thresh, ['Absolute Distance', 'Top Values', 'Top Values (No Thresh)'])
                # plt.savefig('venn_diagram_false_no_thresh.png')
                # plt.close()

                # # Create a Venn Diagram showing overlap between the methods
                # distance_set = set(self.distance_values.keys())
                # top_set = set(self.top_values.keys())
                # rel_dist_set = set(self.relative_distance_values.keys())
                
                # plt.switch_backend('Agg')

                # venn = venn3([distance_set, top_set, rel_dist_set], ('Greatest Distance', 'Top Values', 'Relative Distance'))
                # plt.title('Venn Diagram of Prediction Overlap')

                # # Save the plot to a file
                # plt.savefig('venn_diagram.png')

                # total_known = 0
                # total_unknown = 0
                # total = 0
                # true = 0
                # total_false = 0

                # for k,v in self.distance_values.items():
                #     fam_id = _extract_id(k)
                #     if fam_id in annotation_dict:
                #         if annotation_dict[fam_id] == v[0]:
                #             true +=1
                #         else:
                #             total_false += 1
                #         total_known += 1
                #     else:
                #         total_unknown += 1
                #     total +=1
                
                # print(f"Total Items: {len(self.distance_values)}")
                # print(f"True for Distance: {true}")
                # print(f"False for Distance: {total_false}")
                # print(f"Unknown for Distance: {total_unknown}")
                # print(f"Total for Distance: {total}")
                # print(f"Total_known for top_values_no_thresh: {total_known}")

                # total_known_2 = 0
                # total_2 = 0
                # true_2 = 0
                # total_false_2 = 0
                # total_unknown_2 = 0
                # for k,v in self.top_values.items():
                #     fam_id = _extract_id(k)
                #     if fam_id in annotation_dict:
                #         if annotation_dict[fam_id] == v[0]:
                #             true_2 +=1
                #         else:
                #             total_false_2 += 1
                #         total_known_2 += 1
                #     else:
                #         total_unknown_2 += 1
                #     total_2 +=1

                # print(f"Total Items: {len(self.top_values)}")
                # print(f"True for Top: {true_2}")
                # print(f"False for Top: {total_false_2}")
                # print(f"Unknown for Top: {total_unknown_2}")
                # print(f"Total for Top: {total_2}")
                # print(f"Total_known for Top: {total_known_2}")

                # total_known_3 = 0
                # total_3  = 0
                # true_3 = 0
                # total_false_3 = 0
                # total_unknown_3 = 0
                # for k,v in self.top_values_no_thresh.items():
                #     fam_id = _extract_id(k)
                #     if fam_id in annotation_dict:
                #         if annotation_dict[fam_id] == v[0]:
                #             true_3 += 1
                #         else:
                #             total_false_3 += 1
                #         total_known_3 += 1
                #     else:
                #         total_unknown_3 += 1
                #     total_3 +=1

                # print(f"Total Items: {len(self.top_values_no_thresh)}")
                # print(f"True for top_values_no_thresh: {true_3}")
                # print(f"False for top_values_no_thresh: {total_false_3}")
                # print(f"Unknown for top_values_no_thresh: {total_unknown_3}")
                # print(f"Total for top_values_no_thresh: {total_3}")
                # print(f"Total_known for top_values_no_thresh: {total_known_3}")

                # total_known_4 = 0
                # total_4  = 0
                # true_4 = 0
                # total_false_4 = 0
                # total_unknown_4 = 0
                # for k,v in self.top_values_no_thresh.items():
                #     fam_id = _extract_id(k)
                #     if fam_id in annotation_dict:
                #         if annotation_dict[fam_id] == v[0]:
                #             true_4 += 1
                #         else:
                #             total_false_4 += 1
                #         total_known_4 += 1
                #     else:
                #         total_unknown_4 += 1
                #     total_4 +=1

                # print(f"Total Items: {len(self.relative_distance_values)}")
                # print(f"True for top_values_no_thresh: {true_4}")
                # print(f"False for top_values_no_thresh: {total_false_4}")
                # print(f"Unknown for top_values_no_thresh: {total_unknown_4}")
                # print(f"Total for top_values_no_thresh: {total_4}")
                # print(f"Total_known for top_values_no_thresh: {total_known_4}")


            def execute_all(self):
                """
        Execute the entire sequence of operations in the KmerCompare process.
        """
                self.load_data()
                self.generate_kmer_counts()
                self.construct_kmer_counts_dataframe()
                self.match_kmer_counts_format()
                self.cosine_similarity()
                # self.format_and_write_output()
                self.select_top_no_threshold()
                self.select_top_above_threshold()
                self.select_by_greatest_distance()
                self.select_by_greatest_relative_distance()
                self.calculate_accuracy_greatest_distance()



        apply = KmerCompare(
            input.compare_associations,
            input.data,
            # input.confidence_associations,
            input.decoy_stats,
            input.annotation,
            output.seq_ann,
            output.kmer_summary,
        )
        apply.execute_all()
        skm.utils.log_runtime(log[0], start_time)
