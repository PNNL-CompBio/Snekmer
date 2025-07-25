# ---------------------------------------------------------
# Imports
# ---------------------------------------------------------

import numpy as np
from typing import Tuple, List, Optional, Dict, Any, Sequence, Mapping
import pandas as pd
import itertools
import snekmer as skm

from snekmer.learn import apply_selection_method


# ---------------------------------------------------------
# Files and Parameters
# ---------------------------------------------------------

config = snakemake.config

# ---------------------------------------------------------
# Run script
# ---------------------------------------------------------


class Evaluator:
    """
    The Evaluator class processes predictions and generates confidence metrics.

    Attributes:
    input_data (list): List of paths to input files.
    output_glob (str): Path to save the output global crosstab.
    reverse_decoy_stats (str): Path to family summary statistics for thresholds.
    modifier (float): Weight modifier for confidence merging.
    confidence_data (list, optional): List of paths to base confidence files. Defaults to None.
    """

    input_data: list
    output_glob: str
    reverse_decoy_stats: str
    modifier: float
    confidence_data: None
    true_running_crosstab = None
    false_running_crosstab = None

    def __init__(
        self,
        input_data: list,
        output_glob_path: str,
        reverse_decoy_stats: str,
        modifier: float,
        confidence_data: None,
    ) -> None:
        self.input_data = input_data
        self.output_glob = output_glob_path
        self.reverse_decoy_stats = reverse_decoy_stats
        self.modifier = modifier
        self.confidence_data = confidence_data
        self.true_running_crosstab = None
        self.false_running_crosstab = None

    def read_and_transform_input_data(
        self, file_path: str
    ) -> Tuple[
        List[Tuple[str, str]],
        List[Optional[str]],
        List[float],
        List[str],
        List[str],
        List[str],
        List[float],
    ]:
        """
        Reads and transforms input data from a given CSV file, applying thresholds if necessary.

        Args:
            file_path (str): Path to the input CSV file.

        Returns:
            Tuple containing:
            - top_two (List[Tuple[str, str]]): the top two candidate labels per sequence
            - predictions (List[Optional[str]]): the selected label per sequence
            - deltas (List[float]): the delta scores used in selection
            - result (List[str]): original sequence IDs
            - tf (List[str]): "T" if prediction matches sequence ID, else "F"
            - known (List[str]): "Known" or "Unknown" status per sequence
            - top_scores (List[float]): the score associated with each prediction
        """
        seq_ann_scores = pd.read_csv(
            file_path,
            index_col="__index_level_0__",
            header=0,
            engine="c",
        )
        thresholds_dataframe = pd.read_csv(
            self.reverse_decoy_stats,
            header=0,
            engine="c",
        )
        threshold_type = config["learn_apply"]["threshold"]
        selection_method = config["learn_apply"]["selection"]
        weight_top = config["learn_apply"].get("weight_top", 0.5)
        weight_distance = config["learn_apply"].get("weight_distance", 0.5)

        if threshold_type == "None":
            threshold_type = None
        if threshold_type is not None:
            threshold_dict = dict(
                zip(
                    thresholds_dataframe.family.astype(str),
                    thresholds_dataframe[threshold_type],
                )
            )
        else:
            threshold_dict = None

            self.threshold_dict = threshold_dict

        predictions, deltas, top_two = apply_selection_method(
            seq_ann_scores,
            selection_method,
            threshold_type,
            threshold_dict,
            weight_top,
            weight_distance,
        )

        result = seq_ann_scores.index.tolist()

        tf = []
        known = []
        for i, pred in enumerate(predictions):
            actual = result[i]
            if isinstance(pred, str) and pred in actual:
                tf.append("T")
            else:
                tf.append("F")
            if "unknown" in actual:
                known.append("Unknown")
            else:
                known.append("Known")

        # Retrieve top scores based on predictions
        top_scores = []
        for idx, pred in zip(seq_ann_scores.index, predictions):
            if pred is not None and isinstance(pred, str):
                top_score = seq_ann_scores.at[idx, pred]
            else:
                top_score = np.nan
            top_scores.append(top_score)

        return top_two, predictions, deltas, result, tf, known

    def create_true_false_crosstabs(
        self, diff_dataframe: pd.DataFrame
    ) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """
        Creates crosstabs for True and False predictions based on the given dataframe.

        Args:
            diff_dataframe (DataFrame): DataFrame with differences and metrics.

        Returns:
            tuple: Crosstabs for True and False predictions.
        """
        # Filter out rows where Prediction is None
        valid_predictions = diff_dataframe.dropna(subset=["Prediction"])

        known_true_diff_df = valid_predictions[
            (valid_predictions["Known/Unknown"] == "Known")
            & (valid_predictions["T/F"] == "T")
        ]
        known_false_diff_dataframe = valid_predictions[
            (valid_predictions["Known/Unknown"] == "Known")
            & (valid_predictions["T/F"] == "F")
        ]

        true_crosstab = pd.crosstab(
            known_true_diff_df.Prediction,
            known_true_diff_df.Difference,
        )
        false_crosstab = pd.crosstab(
            known_false_diff_dataframe.Prediction,
            known_false_diff_dataframe.Difference,
        )
        return true_crosstab, false_crosstab

    def handle_running_crosstabs(
        self, true_crosstab: pd.DataFrame, false_crosstab: pd.DataFrame, iteration: int
    ) -> None:
        """
        Handles and updates running crosstabs over iterations.

        Args:
            true_crosstab (DataFrame): Crosstab for True predictions.
            false_crosstab (DataFrame): Crosstab for False predictions.
            iteration (int): Current iteration.
        """
        if iteration == 0:
            self.true_running_crosstab = true_crosstab
            self.false_running_crosstab = false_crosstab
        else:
            self.true_running_crosstab = (
                pd.concat([self.true_running_crosstab, true_crosstab])
                .groupby("Prediction", sort=False)
                .sum(min_count=1)
            ).fillna(0)
            self.false_running_crosstab = (
                pd.concat([self.false_running_crosstab, false_crosstab])
                .groupby("Prediction", sort=False)
                .sum(min_count=1)
            ).fillna(0)

        self.true_running_crosstab, self.false_running_crosstab = (
            self.true_running_crosstab.align(
                self.false_running_crosstab, join="outer", axis=1, fill_value=0
            )
        )
        self.true_running_crosstab.fillna(0, inplace=True)
        self.false_running_crosstab.fillna(0, inplace=True)

    def generate_inputs(self) -> None:
        """
        Format the data so it is ready to be compared to the kmer association matrix.
        """
        for j, f in enumerate(self.input_data):
            (
                two_key_vals,
                predictions,
                deltas,
                result,
                tf,
                known,
            ) = self.read_and_transform_input_data(f)

            diff_dataframe = self.generate_diff_dataframe(
                two_key_vals, predictions, deltas, result, tf, known
            )
            true_crosstab, false_crosstab = self.create_true_false_crosstabs(
                diff_dataframe
            )
            self.handle_running_crosstabs(true_crosstab, false_crosstab, j)

    def generate_global_crosstab(self) -> pd.DataFrame:
        """
        Generates a global crosstab based on the running True and False crosstabs.

        Returns:
            DataFrame: Computed global crosstab.
        """
        return self.true_running_crosstab / (
            self.true_running_crosstab + self.false_running_crosstab
        )

    def calculate_distributions(self) -> Tuple[pd.Series, pd.Series]:
        """
        Calculates distributions for True and False predictions.

        Returns:
            tuple: Total distributions for True and False predictions.
        """
        true_total_dist = self.true_running_crosstab.sum(numeric_only=True, axis=0)
        false_total_dist = self.false_running_crosstab.sum(numeric_only=True, axis=0)

        return true_total_dist, false_total_dist

    def compute_ratio_distribution(
        self, true_total_dist: pd.Series, false_total_dist: pd.Series
    ) -> Tuple[pd.Series, pd.Series, pd.Series]:
        """
        Computes ratio_dist, total_sum, and inter_sum.

        Args:
            true_total_dist (pd.Series): Counts for True predictions at each increment.
            false_total_dist (pd.Series): Counts for False predictions at each increment.

        Returns:
            Tuple[pd.Series, pd.Series, pd.Series]:
                - ratio_dist: ratio of True/(True+False), interpolated over 0.00 to 1.00 by 0.01.
                - total_sum: raw cumulative sum of (True + False) at each increment (missing filled with 0).
                - inter_sum: interpolated version of total_sum over the uniform increments.
        """
        ratio_total_dist = true_total_dist / (true_total_dist + false_total_dist)
        raw_total_sum = (true_total_dist + false_total_dist).copy()
        new_index = pd.Index(
            [round(i, 2) for i in np.arange(0, 1.01, 0.01)],
            name="Difference",
        )
        ratio_total_dist = ratio_total_dist.reindex(new_index)
        ratio_total_dist = ratio_total_dist.interpolate(method="linear")
        ratio_total_dist.fillna(0, inplace=True)
        total_sum = raw_total_sum.reindex(new_index, fill_value=0)
        inter_sum = total_sum.interpolate(method="linear")
        inter_sum.fillna(0, inplace=True)

        return ratio_total_dist, total_sum, inter_sum

    def check_confidence_merge(
        self, new_ratio_dist: pd.Series, total_sum: pd.Series, inter_sum: pd.Series
    ) -> pd.DataFrame:
        """
        Merge the computed distributions with a base confidence file if available.

        Args:
            new_ratio_dist (pd.Series): ratio of True/(True+False) over the uniform increments.
            total_sum (pd.Series): raw cumulative counts of (True + False) at each increment.
            inter_sum (pd.Series): interpolated cumulative counts over the uniform increments.

        Returns:
            pd.DataFrame: DataFrame containing:
                - confidence: merged or new confidence scores per increment
                - weight: combined weight of prior and current data
                - totalSum: updated raw cumulative counts
                - interSum: updated interpolated cumulative counts
                - cur_sum: total sample count across both True and False
        """
        # current_weight is total number of sequences processed
        current_weight = (
            self.true_running_crosstab.values.sum()
            + self.false_running_crosstab.values.sum()
        )

        # Prepare the updated_data DataFrame with the new columns
        updated_data = pd.DataFrame(new_ratio_dist, columns=["confidence"])
        updated_data["weight"] = current_weight
        updated_data["totalSum"] = total_sum  # no interpolation, raw cumulative counts
        updated_data["interSum"] = inter_sum  # interpolated cumulative counts

        if self.confidence_data and len(self.confidence_data) == 1:
            prior_conf = pd.read_csv(self.confidence_data[0], index_col="Difference")
            print(f"Prior Confidence Data:\n{prior_conf}")

            total_weight_prior = prior_conf["weight"]
            k_factor = 1 + self.modifier * (
                current_weight / (current_weight + total_weight_prior)
            )
            out_weight = total_weight_prior + current_weight
            weighted_current = k_factor * current_weight
            total_weight = total_weight_prior + weighted_current
            prior_weighted_score = prior_conf["confidence"] * prior_conf["weight"]
            current_weighted_score = updated_data["confidence"] * weighted_current

            updated_confidence = (
                prior_weighted_score + current_weighted_score
            ) / total_weight
            updated_data["confidence"] = updated_confidence

            # Merge totalSum and interSum with prior if they exist:
            if "totalSum" in prior_conf.columns:
                updated_data["totalSum"] += prior_conf["totalSum"]
            if "interSum" in prior_conf.columns:
                updated_data["interSum"] += prior_conf["interSum"]

                # Update weights and sums if "cur_sum" column is used
            if "cur_sum" in prior_conf.columns:
                # sum might represent cumulative distributions at increments.
                # We can merge them similarly by adding or using another merging strategy:
                updated_data["cur_sum"] = prior_conf["cur_sum"] + (
                    self.true_running_crosstab.sum() + self.false_running_crosstab.sum()
                )
            else:
                # If no prior sum, just create it now
                updated_data["cur_sum"] = (
                    self.true_running_crosstab.sum() + self.false_running_crosstab.sum()
                )

            updated_data["weight"] = out_weight
            print(f"Final Confidence Data\n{updated_data}")
        else:
            print(
                "Base confidence file not found or multiple files present. Only one file is allowed in baseConfidence."
            )
            # If no merging, assign sum as needed
            updated_data["cur_sum"] = (
                self.true_running_crosstab.sum() + self.false_running_crosstab.sum()
            )

        updated_data.fillna(0, inplace=True)
        return updated_data

    def generate_diff_dataframe(
        self,
        two_key_vals: Mapping[str, Sequence[Any]],
        predictions: Sequence[Optional[str]],
        deltas: Sequence[Optional[float]],
        result: Sequence[str],
        tf: Sequence[str],
        known: Sequence[str],
    ) -> pd.DataFrame:
        """
        Generates a dataframe showing the differences and other metrics.

        Args:
            two_key_vals (Mapping[str, Sequence[Any]]):
                Two-key stored values (e.g., top two scores or top score and threshold),
                with keys "key_value_one" and "key_value_two".
            predictions (Sequence[Optional[str]]): Predictions made by the selection method.
            deltas (Sequence[Optional[float]]): Delta values calculated by the selection method.
            result (Sequence[str]): Actual labels.
            tf (Sequence[str]): List of "T"/"F" labels.
            known (Sequence[str]): List of "Known"/"Unknown" labels.

        Returns:
            pd.DataFrame: Constructed dataframe with columns:
                - Top: first key value
                - Second: second key value
                - Difference: rounded delta values
                - Prediction: predicted label
                - Actual: actual label
                - T/F: correctness flag
                - Known/Unknown: knownness flag
        """
        rounded_deltas = [round(num, 2) if num is not None else 0 for num in deltas]

        diff_dataframe = pd.DataFrame(
            {
                "Top": two_key_vals["key_value_one"],
                "Second": two_key_vals["key_value_two"],
                "Difference": rounded_deltas,
                "Prediction": predictions,
                "Actual": result,
                "T/F": tf,
                "Known/Unknown": known,
            }
        )

        tmp = pd.to_numeric(diff_dataframe["Difference"], errors="coerce")
        diff_dataframe["Difference"] = tmp

        return diff_dataframe

    def execute_all(self) -> None:
        """
        Executes all functions in order.

        This method does the following:

        1. Processes each input file by reading, transforming, and storing intermediate variables.
        2. After processing all files, generates a global crosstab to summarize the data.
        3. Calculates distributions for True and False predictions. If base confidence files are provided,
        it merges the calculated distributions with the base ones.
        4. Computes a ratio distribution based on the True and False total distributions.
        5. Finally, writes the ratio distribution and global crosstab to csv.
        """
        self.generate_inputs()
        self.generate_global_crosstab()
        true_dist, false_dist = self.calculate_distributions()
        ratio_dist, total_sum, inter_sum = self.compute_ratio_distribution(
            true_dist, false_dist
        )
        global_confidence = self.check_confidence_merge(
            ratio_dist, total_sum, inter_sum
        )
        global_confidence.to_csv(self.output_glob)


evaluator = Evaluator(
    snakemake.input.eval_apply_data,
    snakemake.output.eval_glob,
    snakemake.input.reverse_decoy_stats,
    snakemake.params.modifier,
    snakemake.input.base_confidence,
)
evaluator.execute_all()
