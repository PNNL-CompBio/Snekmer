# ---------------------------------------------------------
# Imports
# ---------------------------------------------------------

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
# ---------------------------------------------------------
# Files and Parameters
# ---------------------------------------------------------

config = snakemake.config

# change matplotlib backend to non-interactive
plt.switch_backend("Agg")

# ---------------------------------------------------------
# Run script
# ---------------------------------------------------------

class Evaluator:
    """
    The Evaluator class processes predictions and generates confidence metrics.

    Attributes:
    inputData (list): List of paths to input files.
    outputGlobPath (str): Path to save the output global crosstab.
    reverseDecoyStats (str): Path to family summary statistics for thresholds.
    modifier (float): Weight modifier for confidence merging.
    confidenceData (list, optional): List of paths to base confidence files. Defaults to None.
    """

    def __init__(
        self,
        inputData,
        outputGlobPath,
        reverseDecoyStats,
        modifier,
        confidenceData=None,
    ):
        self.inputData = inputData
        self.outputGlob = outputGlobPath
        self.reverseDecoyStats = reverseDecoyStats
        self.modifier = modifier
        self.confidenceData = confidenceData
        self.trueRunningCrosstab = None
        self.falseRunningCrosstab = None

    def readAndTransformInputData(self, file_path):
        """
Reads and transforms input data from a given CSV file, applying thresholds if necessary.

Args:
    file_path (str): Path to the input CSV file.

Returns:
    tuple: Parsed and transformed values from the input data.
"""
        seqAnnScores = pd.read_csv(
            file_path,
            index_col="__index_level_0__",
            header=0,
            engine="c",
        )

        # Load thresholds

        thresholds_df = pd.read_csv(
            self.reverseDecoyStats,
            header=0,
            engine="c",
        )

        threshold_type = config["learnapp"]["threshold"]
        selectionMethod = config["learnapp"]["selection"]

        if threshold_type == "None":
            threshold_type = None
        if threshold_type is not None:
            threshold_dict = dict(
                zip(
                    thresholds_df.Family.astype(str),
                    thresholds_df[threshold_type],
                )
            )

            self.threshold_dict = threshold_dict  # Store for later use

            # Apply the selection method
        predictions, deltas, topTwo = self.apply_selectionMethod(
            seqAnnScores, selectionMethod, threshold_type
        )

        result = seqAnnScores.index.tolist()

        # Generate True/False labels based on predictions
        tf = []
        known = []
        for i, pred in enumerate(predictions):
            actual = result[i]
            if isinstance(pred, str) and pred in actual:
                tf.append("T")
            else:
                tf.append("F")
                # Determine if the sequence is known or unknown
            if "unknown" in actual:
                known.append("Unknown")
            else:
                known.append("Known")

                # Retrieve top scores based on predictions
        top_scores = []
        for idx, pred in zip(seqAnnScores.index, predictions):
            if pred is not None and isinstance(pred, str):
                top_score = seqAnnScores.at[idx, pred]
            else:
                top_score = np.nan
            top_scores.append(top_score)

        return topTwo, predictions, deltas, result, tf, known

    def createTrueFalseCrosstabs(self, diffDataframe):
        """
Creates crosstabs for True and False predictions based on the given dataframe.

Args:
    diffDataframe (DataFrame): DataFrame with differences and metrics.

Returns:
    tuple: Crosstabs for True and False predictions.
"""
        # Filter out rows where Prediction is None
        valid_predictions = diffDataframe.dropna(subset=["Prediction"])

        known_true_diffDataframe = valid_predictions[
            (valid_predictions["Known/Unknown"] == "Known")
            & (valid_predictions["T/F"] == "T")
        ]
        known_false_diffDataframe = valid_predictions[
            (valid_predictions["Known/Unknown"] == "Known")
            & (valid_predictions["T/F"] == "F")
        ]

        trueCrosstab = pd.crosstab(
            known_true_diffDataframe.Prediction,
            known_true_diffDataframe.Difference,
        )
        falseCrosstab = pd.crosstab(
            known_false_diffDataframe.Prediction,
            known_false_diffDataframe.Difference,
        )
        return trueCrosstab, falseCrosstab

    def handleRunningCrosstabs(self, trueCrosstab, falseCrosstab, iteration):
        """
Handles and updates running crosstabs over iterations.

Args:
    trueCrosstab (DataFrame): Crosstab for True predictions.
    falseCrosstab (DataFrame): Crosstab for False predictions.
    iteration (int): Current iteration.
"""
        if iteration == 0:
            self.trueRunningCrosstab = trueCrosstab
            self.falseRunningCrosstab = falseCrosstab
        else:
            self.trueRunningCrosstab = (
                pd.concat([self.trueRunningCrosstab, trueCrosstab])
                .groupby("Prediction", sort=False)
                .sum(min_count=1)
            ).fillna(0)
            self.falseRunningCrosstab = (
                pd.concat([self.falseRunningCrosstab, falseCrosstab])
                .groupby("Prediction", sort=False)
                .sum(min_count=1)
            ).fillna(0)

            # Ensure that both crosstabs have the same columns and index
        self.trueRunningCrosstab, self.falseRunningCrosstab = (
            self.trueRunningCrosstab.align(
                self.falseRunningCrosstab, join="outer", axis=1, fill_value=0
            )
        )
        self.trueRunningCrosstab.fillna(0, inplace=True)
        self.falseRunningCrosstab.fillna(0, inplace=True)

    def generateInputs(self):
        for j, f in enumerate(self.inputData):
            (
                two_key_vals,
                predictions,
                deltas,
                result,
                tf,
                known,
            ) = self.readAndTransformInputData(f)

            diffDataframe = self.generateDiffDataframe(
                two_key_vals, predictions, deltas, result, tf, known
            )
            trueCrosstab, falseCrosstab = self.createTrueFalseCrosstabs(
                diffDataframe
            )
            self.handleRunningCrosstabs(trueCrosstab, falseCrosstab, j)
        return

    def generateGlobalCrosstab(self):
        """
Generates a global crosstab based on the running True and False crosstabs.

Returns:
    DataFrame: Computed global crosstab.
"""
        return self.trueRunningCrosstab / (
            self.trueRunningCrosstab + self.falseRunningCrosstab
        )

    def calculateDistributions(self):
        """
Calculates distributions for True and False predictions.

Returns:
    tuple: Total distributions for True and False predictions.
"""
        trueTotalDist = self.trueRunningCrosstab.sum(numeric_only=True, axis=0)
        false_total_dist = self.falseRunningCrosstab.sum(
            numeric_only=True, axis=0
        )

        return trueTotalDist, false_total_dist

    def computeRatioDistribution(self, trueTotalDist, false_total_dist):
        """
Computes ratioDist, total_sum, and inter_sum.

- ratioDist: as before, a ratio of True/(True+False), interpolated over 0.0 to 1.0 increments.
- total_sum: the raw cumulative sum of (True + False) counts at each increment, without interpolation.
            For increments not originally present, we set 0 to indicate no data (no interpolation).
- inter_sum: an interpolated version of total_sum over the same increments, filling in gaps smoothly.
"""
        # Calculate ratio (as before)
        ratioTotalDist = trueTotalDist / (trueTotalDist + false_total_dist)

        # This is the raw total cumulative sum of samples at each increment
        raw_total_sum = (trueTotalDist + false_total_dist).copy()

        # Define a uniform index from 0.0 to 1.0 with 0.01 increments
        newIndex = pd.Index(
            [round(i, 2) for i in pd.np.arange(0, 1.01, 0.01)],
            name="Difference",
        )

        # Reindex ratioDist to the uniform index, then interpolate to fill missing increments
        ratioTotalDist = ratioTotalDist.reindex(newIndex)
        ratioTotalDist = ratioTotalDist.interpolate(method="linear")
        ratioTotalDist.fillna(0, inplace=True)

        # For total_sum, reindex to the same uniform scale but do not interpolate:
        # Instead of interpolating, we simply fill missing increments with 0.
        total_sum = raw_total_sum.reindex(newIndex, fill_value=0)

        # Now create inter_sum by interpolating total_sum
        inter_sum = total_sum.interpolate(method="linear")
        inter_sum.fillna(0, inplace=True)

        return ratioTotalDist, total_sum, inter_sum

    def checkConfidenceMerge(self, newRatioDist, total_sum, inter_sum):
        """
Merge the computed distributions with a base confidence file if available.
Now we have total_sum and inter_sum to include in the output DataFrame.
"""
        # currentWeight is total number of sequences processed
        currentWeight = (
            self.trueRunningCrosstab.values.sum()
            + self.falseRunningCrosstab.values.sum()
        )

        # Prepare the updatedData DataFrame with the new columns
        updatedData = pd.DataFrame(newRatioDist, columns=["confidence"])
        updatedData["weight"] = currentWeight
        updatedData["total_sum"] = (
            total_sum  # no interpolation, raw cumulative counts
        )
        updatedData["inter_sum"] = inter_sum  # interpolated cumulative counts

        if self.confidenceData and len(self.confidenceData) == 1:
            priorConf = pd.read_csv(
                self.confidenceData[0], index_col="Difference"
            )
            print(f"Prior Confidence Data:\n{priorConf}")

            totalWeightPrior = priorConf["weight"]
            kFactor = 1 + self.modifier * (
                currentWeight / (currentWeight + totalWeightPrior)
            )
            outWeight = totalWeightPrior + currentWeight
            weightedCurrent = kFactor * currentWeight
            totalWeight = totalWeightPrior + weightedCurrent
            priorWeightedScore = priorConf["confidence"] * priorConf["weight"]
            currentWeighted_score = updatedData["confidence"] * weightedCurrent

            updatedConfidence = (
                priorWeightedScore + currentWeighted_score
            ) / totalWeight
            updatedData["confidence"] = updatedConfidence

            # Merge total_sum and inter_sum with prior if they exist:
            if "total_sum" in priorConf.columns:
                updatedData["total_sum"] += priorConf["total_sum"]
            if "inter_sum" in priorConf.columns:
                updatedData["inter_sum"] += priorConf["inter_sum"]

                # Update weights and sums if "cur_sum" column is used
            if "cur_sum" in priorConf.columns:
                # sum might represent cumulative distributions at increments.
                # We can merge them similarly by adding or using another merging strategy:
                updatedData["cur_sum"] = priorConf["cur_sum"] + (
                    self.trueRunningCrosstab.sum()
                    + self.falseRunningCrosstab.sum()
                )
            else:
                # If no prior sum, just create it now
                updatedData["cur_sum"] = (
                    self.trueRunningCrosstab.sum()
                    + self.falseRunningCrosstab.sum()
                )

            updatedData["weight"] = outWeight
            print(f"Final Confidence Data\n{updatedData}")
        else:
            print(
                "Base confidence file not found or multiple files present. Only one file is allowed in baseConfidence."
            )
            # If no merging, assign sum as needed
            updatedData["cur_sum"] = (
                self.trueRunningCrosstab.sum() + self.falseRunningCrosstab.sum()
            )

        updatedData.fillna(0, inplace=True)
        return updatedData

    def saveResults(
        self,
        ratioTotalDist,
        outputPath,
    ):
        """
Saves the computed results to specified paths.

Args:
    ratioTotalDist (Series): Ratio distribution to be saved.
    outputPath (str): Path to save the ratio distribution.

"""
        ratioTotalDist.to_csv(outputPath)
        return

    def executeAll(self):
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
        self.generateInputs()
        ratioCrosstab = self.generateGlobalCrosstab()
        trueDist, falseDist = self.calculateDistributions()
        ratioDist, total_sum, inter_sum = self.computeRatioDistribution(
            trueDist, falseDist
        )
        globalConfidence = self.checkConfidenceMerge(
            ratioDist, total_sum, inter_sum
        )
        self.saveResults(globalConfidence, self.outputGlob)

    def generateDiffDataframe(
        self, two_key_vals, predictions, deltas, result, tf, known
    ):
        """
Generates a dataframe showing the differences and other metrics.

Args:
    two_key_vals (array-like): Two Key stored values, this can be top two scores, or top score and threshold, etc
    predictions (Series): Predictions made by the selection method.
    deltas (Series): Delta values calculated by the selection method.
    result (list): Actual labels.
    tf (list): List of "True/False" labels.
    known (list): List of "Known/Unknown" labels.

Returns:
    DataFrame: Constructed dataframe with computed differences and other metrics.
"""
        # rounded_deltas = [ round(num, 2) for num in deltas ]
        rounded_deltas = [
            round(num, 2) if num is not None else 0 for num in deltas
        ]

        # print(f"two_key_vals in generateDiffDataframe:\n {two_key_vals}")
        diffDataframe = pd.DataFrame(
            {
                "Top": two_key_vals["keyValueOne"],
                "Second": two_key_vals["keyValueTwo"],
                "Difference": rounded_deltas,
                "Prediction": predictions,
                "Actual": result,
                "T/F": tf,
                "Known/Unknown": known,
            }
        )

        tmp = pd.to_numeric(diffDataframe["Difference"], errors="coerce")
        diffDataframe["Difference"] = tmp

        return diffDataframe

    def apply_selectionMethod(
        self, seqAnnScores, selectionMethod, threshold_type
    ):

        if selectionMethod == "top_hit" and threshold_type is None:

            def get_topTwo(row):
                top2 = row.nlargest(2)
                return pd.Series(
                    {
                        "keyValueOne": (
                            top2.iloc[0] if len(top2) > 0 else np.nan
                        ),
                        "keyValueOneHeader": (
                            top2.index[0] if len(top2) > 0 else np.nan
                        ),
                        "keyValueTwo": (
                            top2.iloc[1] if len(top2) > 1 else np.nan
                        ),
                        "keyValueTwoHeader": (
                            top2.index[1] if len(top2) > 1 else np.nan
                        ),
                    }
                )

            keyValsDf = seqAnnScores.apply(get_topTwo, axis=1)
            predictions = keyValsDf["keyValueOneHeader"].tolist()
            deltas = (
                keyValsDf["keyValueOne"] - keyValsDf["keyValueTwo"]
            ).tolist()
            topTwo = keyValsDf[["keyValueOne", "keyValueTwo"]]

            return predictions, deltas, topTwo

        elif selectionMethod == "top_hit" and threshold_type is not None:

            def get_topTwo(row):
                # Filter values based on threshold, keeping only values above the threshold
                thresholds = row.index.map(
                    self.threshold_dict
                ).to_numpy()  # Map thresholds to each family (column) for this row
                filtered_row = row.where(row >= thresholds, np.nan)

                # Get the top two values and their column headers, handling cases with fewer than two valid entries
                top2 = filtered_row.nlargest(2)
                return pd.Series(
                    {
                        "keyValueOne": (
                            top2.iloc[0] if len(top2) > 0 else np.nan
                        ),
                        "keyValueOneHeader": (
                            top2.index[0] if len(top2) > 0 else np.nan
                        ),
                        "keyValueTwo": (
                            top2.iloc[1] if len(top2) > 1 else np.nan
                        ),
                        "keyValueTwoHeader": (
                            top2.index[1] if len(top2) > 1 else np.nan
                        ),
                    }
                )

                # Apply the function to each row of seqAnnScores

            keyValsDf = seqAnnScores.apply(get_topTwo, axis=1)

            # Extracting the predictions and deltas for output
            predictions = keyValsDf["keyValueOneHeader"].tolist()
            deltas = (
                keyValsDf["keyValueOne"] - keyValsDf["keyValueTwo"]
            ).tolist()
            topTwo = keyValsDf[["keyValueOne", "keyValueTwo"]]
            return predictions, deltas, topTwo

            # The top two values compared here are: Greatest dist from threshold and Second greatest dist from Threshold
        elif selectionMethod == "greatest_distance":
            thresholds = seqAnnScores.columns.to_series().map(
                self.threshold_dict
            )
            distances = seqAnnScores.subtract(thresholds, axis=1)
            filtered_distances = distances.where(distances >= 0, np.nan)
            partitioned_distances = np.partition(
                np.nan_to_num(filtered_distances.values, nan=-np.inf),
                -2,
                axis=1,
            )

            topTwoIdx = np.argpartition(
                np.nan_to_num(filtered_distances.values, nan=-np.inf),
                -2,
                axis=1,
            )[:, -2:]

            sortedTopTwoIdx = np.argsort(
                filtered_distances.values[
                    np.arange(filtered_distances.shape[0])[:, None], topTwoIdx
                ],
                axis=1,
            )[:, ::-1]

            topTwoDistances = np.take_along_axis(
                filtered_distances.values,
                topTwoIdx[
                    np.arange(filtered_distances.shape[0])[:, None],
                    sortedTopTwoIdx,
                ],
                axis=1,
            )

            topTwoScores = np.take_along_axis(
                seqAnnScores.values,
                topTwoIdx[
                    np.arange(filtered_distances.shape[0])[:, None],
                    sortedTopTwoIdx,
                ],
                axis=1,
            )

            topTwoHeaders = np.array(filtered_distances.columns)[
                topTwoIdx[
                    np.arange(filtered_distances.shape[0])[:, None],
                    sortedTopTwoIdx,
                ]
            ]

            flattened_headers = topTwoHeaders.flatten()
            flattened_thresholds = thresholds[flattened_headers].values
            topTwoThresholds = flattened_thresholds.reshape(topTwoHeaders.shape)

            keyValsDf = pd.DataFrame(
                {
                    "keyValueOne": topTwoScores[:, 0],
                    "keyValueOneHeader": topTwoHeaders[:, 0],
                    "keyValueOneDistance": topTwoDistances[:, 0],
                    "keyValueOneThreshold": topTwoThresholds[:, 0],
                    "keyValueTwo": topTwoScores[:, 1],
                    "keyValueTwoHeader": topTwoHeaders[:, 1],
                    "keyValueTwoDistance": topTwoDistances[:, 1],
                    "keyValueTwoThreshold": topTwoThresholds[:, 1],
                }
            )

            deltas = (
                keyValsDf["keyValueOneDistance"]
                - keyValsDf["keyValueTwoDistance"]
            ).tolist()
            predictions = keyValsDf["keyValueOneHeader"].tolist()
            topTwo = keyValsDf[["keyValueOne", "keyValueTwo"]]

            return predictions, deltas, topTwo

        elif selectionMethod == "combined_distance":
            # Method 4: Combined method with dynamic delta calculation
            # Set weights
            weight_top = config["learnapp"].get("weight_top", 0.5)
            weight_distance = config["learnapp"].get("weight_distance", 0.5)
            thresholds = seqAnnScores.columns.to_series().map(
                self.threshold_dict
            )
            distances = seqAnnScores - thresholds
            positiveDistances = distances.where(distances >= 0, np.nan)
            combinedScores = (seqAnnScores * weight_top) + (
                distances * weight_distance
            )
            combinedScores = combinedScores.where(
                positiveDistances.notna(), np.nan
            )
            topCombinedScores = combinedScores.max(axis=1)
            predictions = combinedScores.idxmax(axis=1)
            predictions = predictions.where(~topCombinedScores.isna(), None)

            deltas = []
            topTwoList = []

            ## Method 1/2: Top hit with threshold
            filteredScores = seqAnnScores.where(
                seqAnnScores >= thresholds, np.nan
            )
            top_scores_method1 = filteredScores.max(axis=1)
            top_families_method1 = filteredScores.idxmax(axis=1)

            # Get second highest scores
            temp_scores_method1 = filteredScores.apply(
                lambda row: row[row != row.max()], axis=1
            )
            secondScoresMethod1 = temp_scores_method1.max(axis=1)

            ## Method 3: Greatest distance from threshold
            positiveDistancesMethod3 = distances.where(distances >= 0, np.nan)
            topDistancesMethod3 = positiveDistancesMethod3.max(axis=1)
            topFamiliesMethod3 = positiveDistancesMethod3.idxmax(axis=1)

            # Iterate over each sequence to compute delta and topTwo
            for idx in range(len(predictions)):
                pred_family = predictions.iloc[idx]
                if pred_family is None:
                    deltas.append(None)
                    topTwoList.append([np.nan, np.nan])
                    continue

                    # Get top families from other methods
                topFamilyMethod1 = top_families_method1.iloc[idx]
                topFamilyMethod3 = topFamiliesMethod3.iloc[idx]

                if pred_family == topFamilyMethod1:
                    # Use delta from method 1/2: difference between top two scores
                    delta = (
                        top_scores_method1.iloc[idx]
                        - secondScoresMethod1.iloc[idx]
                    )
                    # For topTwo, use top two scores from method 1
                    topTwo = [
                        top_scores_method1.iloc[idx],
                        secondScoresMethod1.iloc[idx],
                    ]
                elif pred_family == topFamilyMethod3:
                    # Use delta from method 3: difference between top score and threshold
                    original_score = seqAnnScores.loc[
                        seqAnnScores.index[idx], pred_family
                    ]
                    threshold = thresholds[pred_family]
                    delta = original_score - threshold
                    # For topTwo, use original score and threshold
                    topTwo = [original_score, threshold]
                else:
                    # Use delta from combined scores: difference between top two combined scores
                    row_combinedScores = combinedScores.iloc[idx]
                    # Exclude the top prediction to find the second highest combined score
                    secondCombinedScore = row_combinedScores.drop(
                        pred_family
                    ).max()
                    delta = topCombinedScores.iloc[idx] - secondCombinedScore
                    # For topTwo, use top two combined scores
                    topTwo = [topCombinedScores.iloc[idx], secondCombinedScore]

                deltas.append(delta)
                topTwoList.append(topTwo)

                # Create a DataFrame for topTwo
            topTwo = pd.DataFrame(
                topTwoList, columns=["keyValueOne", "keyValueTwo"]
            )

            # Convert predictions to list
            predictions = predictions.tolist()

            return predictions, deltas, topTwo

        else:
            raise ValueError(f"Invalid selection method: {selectionMethod}")


evaluator = Evaluator(
    snakemake.input.evalApplyData,
    snakemake.output.evalGlob,
    snakemake.input.reverseDecoyStats,
    snakemake.params.modifier,
    snakemake.input.baseConfidence,
)
evaluator.executeAll()