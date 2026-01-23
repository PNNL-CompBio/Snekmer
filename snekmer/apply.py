# snekmer/apply.py

# ---------------------------------------------------------
# Imports
# ---------------------------------------------------------

import itertools
import sys

import math
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.csv as csv
from sklearn.metrics.pairwise import cosine_similarity

from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import snekmer as skm
from snekmer.learn import apply_selection_method, generate_kmer_counts, match_kmer_counts_format

# ---------------------------------------------------------
# Classes / Functions
# ---------------------------------------------------------

class KmerCompare:
    """
    Compares k-mer counts across sequences and applies a selection method to generate predictions.

    Args:
        compare_associations (str): Path to the CSV with k-mer totals per sequence.
        data (str): Path to the NPZ file containing raw k-mer data and sequence IDs.
        confidence_associations (str): Path to CSV with global confidence scores.
        decoy_stats (str): Path to CSV with decoy stats for threshold building.
        output_seq_ann (str): Path to write formatted k-mer association output.
        output_kmer_summary (str): Path to write the final prediction summary CSV.
        selection_type (str): Name of the selection method to apply.
        threshold_type (str): Column name in decoy_stats to use for thresholding.
    """
    compare_associations: str
    data: str
    confidence_associations: str
    decoy_stats: str
    output_seq_ann: str
    output_kmer_summary: str
    selection_type: str
    threshold_ype: str
    df: Optional[pd.DataFrame]
    
    
    def __init__(
        self,
        compare_associations: str,
        data: str,
        confidence_associations: str,
        decoy_stats: str,
        output_seq_ann: str,
        output_kmer_summary: str,
        selection_type: str,
        threshold_type: str,
    )  -> None:

        self.compare_associations = compare_associations
        self.data = data
        self.confidence_associations = confidence_associations
        self.decoy_stats = decoy_stats
        self.output_seq_ann = output_seq_ann
        self.output_kmer_summary = output_kmer_summary
        self.selection_type = selection_type
        self.threshold_type = threshold_type
        self.df = None
        self.kmer_totals = None
        self.seq_kmer_dict = {}

    def load_data(self) -> None:
        """
        Loads k-mer totals and sequence data.

        Reads the compare_associations CSV into a pandas DataFrame and
        loads the k-mer list and sequence DataFrame from the NPZ file.
        """
        self.kmer_count_totals = pd.read_csv(
            str(self.compare_associations),
            index_col="__index_level_0__",
            header=0,
            engine="c",
        )
        kmer_list, df = skm.io.load_npz(self.data)
        self.kmer_list = kmer_list[0]
        self.seq_ids = df["sequence_id"]
        self.df = df

    def construct_kmer_counts_dataframe(self) -> None:
        """
        Constructs a DataFrame of k-mer counts across sequences.

        Transforms the per-sequence k-mer counts dictionary into a DataFrame,
        inserts total and per-sequence counts in the appropriate format.
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
        self.kmer_counts.columns = ["Sequence count"] + list(self.kmer_list)
        self.kmer_counts.index = ["Totals"] + list(self.seq_kmer_dict.keys())


    def apply_cosine_similarity(self) -> None:
        """
        Computes cosine similarity between k-mer totals and per-sequence counts using sklearn.

        Updates `kmer_count_totals` to the matrix of pairwise cosine similarities.
        """
        cosine_dataframe = cosine_similarity(
            self.kmer_count_totals, self.kmer_counts
        ).T
        self.kmer_count_totals = pd.DataFrame(
            cosine_dataframe,
            columns=self.kmer_count_totals.index,
            index=self.kmer_counts.index,
        )
    
    def format_and_write_output(self, save_apply_associations: str) -> None:
        """
        Formats and writes the final Apply output files.

        Args:
            save_apply_associations (str): If truthy, writes the k-mer association
                DataFrame to `output_seq_ann`.
        """
        if save_apply_associations:
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
            else:
                prediction, score, delta = None, None, None

            if (score is None) or (isinstance(score, float) and (np.isnan(score) or np.isinf(score))):
                score = 0.0

            if delta is None or (isinstance(delta, float) and (np.isnan(delta) or np.isinf(delta))):
                delta = 0.0

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


    def build_threshold_dict(self) -> None:
        """
        Builds a threshold lookup dictionary from decoy stats.

        Reads the `decoy_stats` CSV and maps each Family to its threshold value
        based on `threshold_type`.

        """
        if self.threshold_type is None:
            self.threshold_dict = {}
        else:
            df = pd.read_csv(str(self.decoy_stats), engine="c")
            self.threshold_dict = dict(zip(df.family.astype(str),
                                        df[self.threshold_type]))
            
                        
    def execute_all(self, weight_top: str, weight_distance: str, save_apply_associations:str) -> None:
        """
        Runs the full Apply pipeline.

        Loads data, generates k-mer counts, matches formats, computes similarities,
        builds thresholds, applies selection, and writes outputs.
        
        Every sequence will have a prediction. If the score is zero, it indicates no signal though.

        Args:
            weight_top (str): Weight parameter for the top-hit component of selection.
            weight_distance (str): Weight parameter for the distance component of selection.
            save_apply_associations (str): If truthy, saves intermediate associations.
        """
        self.load_data()
        
        self.kmer_list, self.kmer_totals, self.seq_kmer_dict = generate_kmer_counts(self.data, 
                                                                            self.kmer_list,
                                                                            self.kmer_totals, 
                                                                            self.seq_kmer_dict, 
                                                                            False)
        self.construct_kmer_counts_dataframe()
        self.kmer_counts, self.kmer_count_totals = match_kmer_counts_format(
            self.kmer_counts, self.kmer_count_totals)
        
        self.kmer_count_totals = self.kmer_count_totals.apply(pd.to_numeric, errors="coerce").fillna(0)
        self.kmer_counts       = self.kmer_counts.apply(pd.to_numeric, errors="coerce").fillna(0)

        self.apply_cosine_similarity()
        
        mat = (
            self.kmer_count_totals
            .replace([np.inf, -np.inf], np.nan)
            .fillna(0)
        )
        # no signal if the best similarity is 0
        self.no_signal_idx = set(mat.index[mat.max(axis=1) <= 0].tolist())

        self.build_threshold_dict()
        
        preds, deltas, top_two = apply_selection_method(
            self.kmer_count_totals,
            self.selection_type,
            self.threshold_type,
            self.threshold_dict,
            weight_top,
            weight_distance)
        
            
        top_two.index = self.kmer_count_totals.index
        self.selected_values = {}
        for seq_id, pred, delta in zip(self.kmer_count_totals.index, preds, deltas):
            if seq_id in self.no_signal_idx:
                self.selected_values[seq_id] = (None, None, None)
                continue

            if pred is None:
                self.selected_values[seq_id] = (None, None, None)
            else:
                top_val = float(top_two.loc[seq_id, "key_value_one"])
                self.selected_values[seq_id] = (pred, top_val, float(0 if delta is None else delta))
                
        self.format_and_write_output(save_apply_associations)
