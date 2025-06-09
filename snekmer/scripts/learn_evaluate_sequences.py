# ---------------------------------------------------------
# Imports
# ---------------------------------------------------------

import numpy as np
from typing import Tuple, List, Optional, Dict, Any, Sequence, Mapping
import pandas as pd
import itertools
import snekmer as skm


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
    inputData (list): List of paths to input files.
    outputGlobPath (str): Path to save the output global crosstab.
    reverseDecoyStats (str): Path to family summary statistics for thresholds.
    modifier (float): Weight modifier for confidence merging.
    confidenceData (list, optional): List of paths to base confidence files. Defaults to None.
    """
    inputData: list
    outputGlob: str
    reverseDecoyStats: str
    modifier: float
    confidenceData: None
    trueRunningCrosstab = None
    falseRunningCrosstab = None

    def __init__(
        self,
        inputData: list,
        outputGlobPath: str,
        reverseDecoyStats: str,
        modifier: float,
        confidenceData: None,
        ) -> None:
        
        self.inputData = inputData
        self.outputGlob = outputGlobPath
        self.reverseDecoyStats = reverseDecoyStats
        self.modifier = modifier
        self.confidenceData = confidenceData
        self.trueRunningCrosstab = None
        self.falseRunningCrosstab = None

    def readAndTransformInputData(self, filePath: str) -> Tuple[
        List[Tuple[str, str]],
        List[Optional[str]],
        List[float],
        List[str],
        List[str],
        List[str],
        List[float]
        ]:
        """
        Reads and transforms input data from a given CSV file, applying thresholds if necessary.

        Args:
            filePath (str): Path to the input CSV file.

        Returns:
            Tuple containing:
            - topTwo (List[Tuple[str, str]]): the top two candidate labels per sequence
            - predictions (List[Optional[str]]): the selected label per sequence
            - deltas (List[float]): the delta scores used in selection
            - result (List[str]): original sequence IDs
            - tf (List[str]): "T" if prediction matches sequence ID, else "F"
            - known (List[str]): "Known" or "Unknown" status per sequence
            - topScores (List[float]): the score associated with each prediction
        """
        
        seqAnnScores = pd.read_csv(
            filePath,
            index_col="__index_level_0__",
            header=0,
            engine="c",
        )

        # Load thresholds
        thresholdsDataframe = pd.read_csv(
            self.reverseDecoyStats,
            header=0,
            engine="c",
        )

        thresholdType = config["learnapp"]["threshold"]
        selectionMethod = config["learnapp"]["selection"]

        if thresholdType == "None":
            thresholdType = None
        if thresholdType is not None:
            thresholdDict = dict(
                zip(
                    thresholdsDataframe.Family.astype(str),
                    thresholdsDataframe[thresholdType],
                )
            )

            self.thresholdDict = thresholdDict  # Store for later use

        # Apply the selection method
        predictions, deltas, topTwo = self.applySelectionMethod(
            seqAnnScores, selectionMethod, thresholdType
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
        topScores = []
        for idx, pred in zip(seqAnnScores.index, predictions):
            if pred is not None and isinstance(pred, str):
                topScore = seqAnnScores.at[idx, pred]
            else:
                topScore = np.nan
            topScores.append(topScore)

        return topTwo, predictions, deltas, result, tf, known

    def createTrueFalseCrosstabs(self, diffDataframe: pd.DataFrame) -> Tuple[
        pd.DataFrame, 
        pd.DataFrame
        ]:
        """
        Creates crosstabs for True and False predictions based on the given dataframe.

        Args:
            diffDataframe (DataFrame): DataFrame with differences and metrics.

        Returns:
            tuple: Crosstabs for True and False predictions.
        """
        # Filter out rows where Prediction is None
        validPredictions = diffDataframe.dropna(subset=["Prediction"])

        knownTrueDiffDataframe = validPredictions[
            (validPredictions["Known/Unknown"] == "Known")
            & (validPredictions["T/F"] == "T")
        ]
        knownFalseDiffDataframe = validPredictions[
            (validPredictions["Known/Unknown"] == "Known")
            & (validPredictions["T/F"] == "F")
        ]

        trueCrosstab = pd.crosstab(
            knownTrueDiffDataframe.Prediction,
            knownTrueDiffDataframe.Difference,
        )
        falseCrosstab = pd.crosstab(
            knownFalseDiffDataframe.Prediction,
            knownFalseDiffDataframe.Difference,
        )
        return trueCrosstab, falseCrosstab

    def handleRunningCrosstabs(self, 
        trueCrosstab: pd.DataFrame, 
        falseCrosstab: pd.DataFrame, 
        iteration: int) -> None:
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

    def generateInputs(self) -> None:
        """
        Format the data so it is ready to be compared to the kmer association matrix.
        """
        for j, f in enumerate(self.inputData):
            (
                twoKeyVals,
                predictions,
                deltas,
                result,
                tf,
                known,
            ) = self.readAndTransformInputData(f)

            diffDataframe = self.generateDiffDataframe(
                twoKeyVals, predictions, deltas, result, tf, known
            )
            trueCrosstab, falseCrosstab = self.createTrueFalseCrosstabs(
                diffDataframe
            )
            self.handleRunningCrosstabs(trueCrosstab, falseCrosstab, j)

    def generateGlobalCrosstab(self) -> pd.DataFrame:
        """
        Generates a global crosstab based on the running True and False crosstabs.

        Returns:
            DataFrame: Computed global crosstab.
        """
        return self.trueRunningCrosstab / (
            self.trueRunningCrosstab + self.falseRunningCrosstab
        )

    def calculateDistributions(self) -> Tuple[pd.Series, pd.Series]:
        """
        Calculates distributions for True and False predictions.

        Returns:
            tuple: Total distributions for True and False predictions.
        """
        trueTotalDist = self.trueRunningCrosstab.sum(numeric_only=True, axis=0)
        falseTotalDist = self.falseRunningCrosstab.sum(
            numeric_only=True, axis=0
        )

        return trueTotalDist, falseTotalDist

    def computeRatioDistribution(
        self,
        trueTotalDist: pd.Series,
        falseTotalDist: pd.Series
        ) -> Tuple[pd.Series, 
               pd.Series, 
               pd.Series]:
        """
        Computes ratioDist, totalSum, and interSum.

        Args:
            trueTotalDist (pd.Series): Counts for True predictions at each increment.
            falseTotalDist (pd.Series): Counts for False predictions at each increment.

        Returns:
            Tuple[pd.Series, pd.Series, pd.Series]:
                - ratioDist: ratio of True/(True+False), interpolated over 0.00 to 1.00 by 0.01.
                - totalSum: raw cumulative sum of (True + False) at each increment (missing filled with 0).
                - interSum: interpolated version of totalSum over the uniform increments.
        """
        ratioTotalDist = trueTotalDist / (trueTotalDist + falseTotalDist)
        # This is the raw total cumulative sum of samples at each increment
        rawTotalSum = (trueTotalDist + falseTotalDist).copy()
        # Define a uniform index from 0.0 to 1.0 with 0.01 increments
        newIndex = pd.Index(
            [round(i, 2) for i in pd.np.arange(0, 1.01, 0.01)],
            name="Difference",
        )
        # Reindex ratioDist to the uniform index, then interpolate to fill missing increments
        ratioTotalDist = ratioTotalDist.reindex(newIndex)
        ratioTotalDist = ratioTotalDist.interpolate(method="linear")
        ratioTotalDist.fillna(0, inplace=True)
        # For totalSum, reindex to the same uniform scale but do not interpolate:
        # Instead of interpolating, we simply fill missing increments with 0.
        totalSum = rawTotalSum.reindex(newIndex, fill_value=0)
        
        # Now create interSum by interpolating totalSum
        interSum = totalSum.interpolate(method="linear")
        interSum.fillna(0, inplace=True)

        return ratioTotalDist, totalSum, interSum

    def checkConfidenceMerge(
        self,
        newRatioDist: pd.Series,
        totalSum: pd.Series,
        interSum: pd.Series
        ) -> pd.DataFrame:
        """
        Merge the computed distributions with a base confidence file if available.

        Args:
            newRatioDist (pd.Series): ratio of True/(True+False) over the uniform increments.
            totalSum (pd.Series): raw cumulative counts of (True + False) at each increment.
            interSum (pd.Series): interpolated cumulative counts over the uniform increments.

        Returns:
            pd.DataFrame: DataFrame containing:
                - confidence: merged or new confidence scores per increment
                - weight: combined weight of prior and current data
                - totalSum: updated raw cumulative counts
                - interSum: updated interpolated cumulative counts
                - cur_sum: total sample count across both True and False
        """
        # currentWeight is total number of sequences processed
        currentWeight = (
            self.trueRunningCrosstab.values.sum()
            + self.falseRunningCrosstab.values.sum()
        )

        # Prepare the updatedData DataFrame with the new columns
        updatedData = pd.DataFrame(newRatioDist, columns=["confidence"])
        updatedData["weight"] = currentWeight
        updatedData["totalSum"] = (
            totalSum  # no interpolation, raw cumulative counts
        )
        updatedData["interSum"] = interSum  # interpolated cumulative counts

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
            currentWeightedScore = updatedData["confidence"] * weightedCurrent

            updatedConfidence = (
                priorWeightedScore + currentWeightedScore
            ) / totalWeight
            updatedData["confidence"] = updatedConfidence

            # Merge totalSum and interSum with prior if they exist:
            if "totalSum" in priorConf.columns:
                updatedData["totalSum"] += priorConf["totalSum"]
            if "interSum" in priorConf.columns:
                updatedData["interSum"] += priorConf["interSum"]

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

    def generateDiffDataframe(
        self,
        twoKeyVals: Mapping[str, Sequence[Any]],
        predictions: Sequence[Optional[str]],
        deltas: Sequence[Optional[float]],
        result: Sequence[str],
        tf: Sequence[str],
        known: Sequence[str]
        ) -> pd.DataFrame:
        """
        Generates a dataframe showing the differences and other metrics.

        Args:
            twoKeyVals (Mapping[str, Sequence[Any]]): 
                Two-key stored values (e.g., top two scores or top score and threshold),
                with keys "keyValueOne" and "keyValueTwo".
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
        roundedDeltas = [
            round(num, 2) if num is not None else 0 for num in deltas
        ]

        diffDataframe = pd.DataFrame(
            {
                "Top": twoKeyVals["keyValueOne"],
                "Second": twoKeyVals["keyValueTwo"],
                "Difference": roundedDeltas,
                "Prediction": predictions,
                "Actual": result,
                "T/F": tf,
                "Known/Unknown": known,
            }
        )

        tmp = pd.to_numeric(diffDataframe["Difference"], errors="coerce")
        diffDataframe["Difference"] = tmp

        return diffDataframe

    def applySelectionMethod(
        self,
        seqAnnScores: pd.DataFrame,
        selectionMethod: str,
        thresholdType: Optional[str]
    ) -> Tuple[List[Optional[str]], List[float], pd.DataFrame]:
        """
        Apply the chosen selection method to sequence annotation scores.

        Args:
            seqAnnScores (pd.DataFrame): DataFrame of annotation scores, indexed by sequence ID.
            selectionMethod (str): Which method to use (e.g., "top_hit").
            thresholdType (Optional[str]): Name of the threshold column, or None.

        Returns:
            Tuple containing:
            - predictions (List[Optional[str]]): One predicted label per sequence.
            - deltas (List[float]): Difference between top two scores per sequence.
            - topTwo (pd.DataFrame): DataFrame with columns "keyValueOne" and "keyValueTwo" of the top two scores.
        """
            
        if selectionMethod == "top_hit" and thresholdType is None:

            def getTopTwo(row: pd.Series) -> pd.Series:
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

            keyValsDf = seqAnnScores.apply(getTopTwo, axis=1)
            predictions = keyValsDf["keyValueOneHeader"].tolist()
            deltas = (
                keyValsDf["keyValueOne"] - keyValsDf["keyValueTwo"]
            ).tolist()
            topTwo = keyValsDf[["keyValueOne", "keyValueTwo"]]

        elif selectionMethod == "top_hit" and thresholdType is not None:

            def getTopTwo(row: pd.Series) -> pd.Series:
                # Filter values based on threshold, keeping only values above the threshold
                thresholds = row.index.map(
                    self.thresholdDict
                ).to_numpy()  # Map thresholds to each family (column) for this row
                filteredRow = row.where(row >= thresholds, np.nan)

                # Get the top two values and their column headers, handling cases with fewer than two valid entries
                top2 = filteredRow.nlargest(2)
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

            keyValsDf = seqAnnScores.apply(getTopTwo, axis=1)

            # Extracting the predictions and deltas for output
            predictions = keyValsDf["keyValueOneHeader"].tolist()
            deltas = (
                keyValsDf["keyValueOne"] - keyValsDf["keyValueTwo"]
            ).tolist()
            topTwo = keyValsDf[["keyValueOne", "keyValueTwo"]]

        # The top two values compared here are: Greatest dist from threshold and Second greatest dist from Threshold
        elif selectionMethod == "greatest_distance":
            thresholds = seqAnnScores.columns.to_series().map(
                self.thresholdDict
            )
            distances = seqAnnScores.subtract(thresholds, axis=1)
            filteredDistances = distances.where(distances >= 0, np.nan)
            partitioned_distances = np.partition(
                np.nan_to_num(filteredDistances.values, nan=-np.inf),
                -2,
                axis=1,
            )

            topTwoIdx = np.argpartition(
                np.nan_to_num(filteredDistances.values, nan=-np.inf),
                -2,
                axis=1,
            )[:, -2:]

            sortedTopTwoIdx = np.argsort(
                filteredDistances.values[
                    np.arange(filteredDistances.shape[0])[:, None], topTwoIdx
                ],
                axis=1,
            )[:, ::-1]

            topTwoDistances = np.take_along_axis(
                filteredDistances.values,
                topTwoIdx[
                    np.arange(filteredDistances.shape[0])[:, None],
                    sortedTopTwoIdx,
                ],
                axis=1,
            )

            topTwoScores = np.take_along_axis(
                seqAnnScores.values,
                topTwoIdx[
                    np.arange(filteredDistances.shape[0])[:, None],
                    sortedTopTwoIdx,
                ],
                axis=1,
            )

            topTwoHeaders = np.array(filteredDistances.columns)[
                topTwoIdx[
                    np.arange(filteredDistances.shape[0])[:, None],
                    sortedTopTwoIdx,
                ]
            ]

            flattenedHeaders = topTwoHeaders.flatten()
            flattenedThresholds = thresholds[flattenedHeaders].values
            topTwoThresholds = flattenedThresholds.reshape(topTwoHeaders.shape)

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

        elif selectionMethod == "combined_distance":
            weightTop = config["learnapp"].get("weight_top", 0.5)
            weightDistance = config["learnapp"].get("weight_distance", 0.5)
            thresholds = seqAnnScores.columns.to_series().map(
                self.thresholdDict
            )
            distances = seqAnnScores - thresholds
            positiveDistances = distances.where(distances >= 0, np.nan)
            combinedScores = (seqAnnScores * weightTop) + (
                distances * weightDistance
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
            topScoresMethodOne = filteredScores.max(axis=1)
            topFamiliesMethodOne = filteredScores.idxmax(axis=1)

            # Get second highest scores
            tempScoresMethodOne = filteredScores.apply(
                lambda row: row[row != row.max()], axis=1
            )
            secondScoresMethodOne = tempScoresMethodOne.max(axis=1)

            ## Method 3: Greatest distance from threshold
            positiveDistancesMethodThree = distances.where(distances >= 0, np.nan)
            topDistancesMethodThree = positiveDistancesMethodThree.max(axis=1)
            topFamiliesMethodThree = positiveDistancesMethodThree.idxmax(axis=1)

            # Iterate over each sequence to compute delta and topTwo
            for idx in range(len(predictions)):
                predFamily = predictions.iloc[idx]
                if predFamily is None:
                    deltas.append(None)
                    topTwoList.append([np.nan, np.nan])
                    continue

                    # Get top families from other methods
                topFamilyMethod1 = topFamiliesMethodOne.iloc[idx]
                topFamilyMethodThree = topFamiliesMethodThree.iloc[idx]

                if predFamily == topFamilyMethod1:
                    delta = (
                        topScoresMethodOne.iloc[idx]
                        - secondScoresMethodOne.iloc[idx]
                    )
                    topTwo = [
                        topScoresMethodOne.iloc[idx],
                        secondScoresMethodOne.iloc[idx],
                    ]
                elif predFamily == topFamilyMethodThree:
                    originalScore = seqAnnScores.loc[
                        seqAnnScores.index[idx], predFamily
                    ]
                    threshold = thresholds[predFamily]
                    delta = originalScore - threshold
                    topTwo = [originalScore, threshold]
                else:
                    rowCombinedScores = combinedScores.iloc[idx]
                    secondCombinedScore = rowCombinedScores.drop(
                        predFamily
                    ).max()
                    delta = topCombinedScores.iloc[idx] - secondCombinedScore
                    topTwo = [topCombinedScores.iloc[idx], secondCombinedScore]

                deltas.append(delta)
                topTwoList.append(topTwo)

                # Create a DataFrame for topTwo
            topTwo = pd.DataFrame(
                topTwoList, columns=["keyValueOne", "keyValueTwo"]
            )

            # Convert predictions to list
            predictions = predictions.tolist()
        else:
            raise ValueError(f"Invalid selection method: {selectionMethod}")
        
        return predictions, deltas, topTwo

    def saveResults(
        self,
        ratioTotalDist: pd.Series,
        outputPath: str
        ) -> None:
        """
        Saves the computed ratio distribution to a CSV file.

        Args:
            ratioTotalDist (pd.Series): Ratio distribution indexed by 'Difference'.
            outputPath (str): File path where the CSV will be written.

        Returns:
            None
        """
        ratioTotalDist.to_csv(outputPath)

    def executeAll(self) -> None:
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
        ratioDist, totalSum, interSum = self.computeRatioDistribution(
            trueDist, falseDist
        )
        globalConfidence = self.checkConfidenceMerge(
            ratioDist, totalSum, interSum
        )
        self.saveResults(globalConfidence, self.outputGlob)
        
evaluator = Evaluator(
    snakemake.input.evalApplyData,
    snakemake.output.evalGlob,
    snakemake.input.reverseDecoyStats,
    snakemake.params.modifier,
    snakemake.input.baseConfidence,
)
evaluator.executeAll()