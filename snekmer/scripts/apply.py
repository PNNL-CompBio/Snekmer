# ---------------------------------------------------------
# Imports
# ---------------------------------------------------------

import sys
import itertools
from typing import Optional
import pandas as pd
import numpy as np
import snekmer as skm
import pyarrow as pa
import pyarrow.csv as csv
from sklearn.metrics.pairwise import cosine_similarity

# ---------------------------------------------------------
# Files and Parameters
# ---------------------------------------------------------
config = snakemake.config


outDir = skm.io.define_output_dir(
    config["alphabet"], config["k"], nested=config["nested_output"]
)

# ---------------------------------------------------------
# Run script
# ---------------------------------------------------------


class KmerCompare:
    """
    Initialize KmerCompare with necessary file paths.

    Args:
        compareAssociations (str): Path to compare associations file.
        data (str): Path to data file.
        confidenceAssociations (str): Path to confidence associations file.
        decoyStats (str): Path to decoy stats file.
        annotation (str): Path to annotation file.
        outputSeqAnn (str): Path to output sequence annotation file.
        outputKmerSummary (str): Path to output kmer summary file.
        selectionType (str): Method selection based on config.
        thresholdType (str): Threshold type from config.
    """
    compareAssociations: str
    data: str
    confidenceAssociations: str
    decoyStats: str
    annotation: str
    outputSeqAnn: str
    outputKmerSummary: str
    selectionType: str
    thresholdType: str
    df: Optional[pd.DataFrame]
    
    def __init__(
        self,
        compareAssociations: str,
        data: str,
        confidenceAssociations: str,
        decoyStats: str,
        annotation: str,
        outputSeqAnn: str,
        outputKmerSummary: str,
        selectionType: str,
        thresholdType: str,
    )  -> None:

        self.compareAssociations = compareAssociations
        self.data = data
        self.confidenceAssociations = confidenceAssociations
        self.annotation = annotation
        self.decoyStats = decoyStats
        self.outputSeqAnn = outputSeqAnn
        self.outputKmerSummary = outputKmerSummary
        self.selectionType = selectionType
        self.thresholdType = thresholdType
        self.df = None

    def loadData(self) -> None:
        """
        Load kmer counts and sequence data from provided files.
        """
        self.kmerCountTotals = pd.read_csv(
            str(self.compareAssociations),
            index_col="__index_level_0__",
            header=0,
            engine="c",
        )
        kmerList, df = skm.io.load_npz(self.data)
        self.kmerList = kmerList[0]
        self.seqIDs = df["sequence_id"]
        self.df = df

    def generateKmerCounts(self) -> None:
        """
        Generate k-mer counts for sequences present in the data.
        """
        self.kmerTotals = [0 for _ in self.kmerList]
        kmerLen = len(self.kmerList[0])
        self.seqKmerDict = {}
        for i, seq in enumerate(self.seqIDs):
            v = self.df["sequence"][i]
            kmerCounts = dict()
            items = [
                v[item : item + kmerLen] for item in range(0, len(v) - kmerLen + 1)
            ]
            for j in items:
                kmerCounts[j] = kmerCounts.get(j, 0) + 1
            store = [kmerCounts.get(item, 0) for item in self.kmerList]
            for i, item in enumerate(self.kmerList):
                self.kmerTotals[i] += kmerCounts.get(item, 0)
            self.seqKmerDict[seq] = store

    def constructKmerCountsDataframe(self) -> None:
        """
        Construct a DataFrame to represent k-mer counts across sequences.
        """
        totalSeqs = len(self.seqKmerDict)
        self.kmerCounts = pd.DataFrame(self.seqKmerDict.values())
        self.kmerCounts.insert(0, "Annotations", 1, True)
        self.kmerTotals.insert(0, totalSeqs)
        self.kmerCounts = pd.DataFrame(
            np.insert(
                self.kmerCounts.values, 0, values=self.kmerTotals, axis=0
            )
        )
        self.kmerCounts.columns = ["Sequence count"] + list(self.kmerList)
        self.kmerCounts.index = ["Totals"] + list(self.seqKmerDict.keys())

    def matchKmerCountsFormat(self) -> None:
        """
        Ensure that the format of the k-mer counts DataFrame matches the expected format.
        """
        if len(str(self.kmerCounts.columns.values[10])) == len(
            str(self.kmerCountTotals.columns.values[10])
        ):
            compareCheck = True
        else:
            compareCheck = False

        if compareCheck:
            check = len(self.kmerCounts.columns.values)
            alphabetInitial = set(
                itertools.chain(
                    *[
                        list(x)
                        for x in self.kmerCounts.columns.values[10:check]
                    ]
                )
            )
            alphabetCompare = set(
                itertools.chain(
                    *[
                        list(x)
                        for x in self.kmerCountTotals.columns.values[
                            10:check
                        ]
                    ]
                )
            )
            if alphabetCompare != alphabetInitial:
                compareCheck = False

        if not compareCheck:
            print("Compare Check Failed. ")
            sys.exit()

        self.kmerCounts.drop("Totals", axis=0, inplace=True)
        self.kmerCounts.drop("Sequence count", axis=1, inplace=True)
        self.kmerCountTotals.drop("Totals", axis=0, inplace=True)
        self.kmerCountTotals.drop("Kmer Count", axis=1, inplace=True)
        self.kmerCountTotals.drop("Sequence count", axis=1, inplace=True)

        columnOrder = list(
            set(self.kmerCounts.columns) | set(self.kmerCountTotals.columns)
        )
        self.kmerCounts = self.kmerCounts.reindex(
            columns=columnOrder, fillValue=0
        )
        self.kmerCountTotals = self.kmerCountTotals.reindex(
            columns=columnOrder, fillValue=0
        )

    def cosineSimilarity(self) -> None:
        """
        Compute cosine similarity between kmer counts of sequences.
        """
        cosineDataframe = cosine_similarity(
            self.kmerCountTotals, self.kmerCounts
        ).T
        self.kmerCountTotals = pd.DataFrame(
            cosineDataframe,
            columns=self.kmerCountTotals.index,
            index=self.kmerCounts.index,
        )

        # Method 0: Hit Hit No Threshold

    def selectTopNoThreshold(self) -> None:
        """
        Select the top hit with no family specific threshold specified
        """
        self.selectedValues = {}
        for rowID, row in self.kmerCountTotals.iterrows():
            if not row.empty:
                sortedRow = row.sort_values(ascending=False)
                topValue = sortedRow.iloc[0]
                topFamily = sortedRow.index[0]
                if len(sortedRow) > 1:
                    secondValue = sortedRow.iloc[1]
                    delta = topValue - secondValue
                else:
                    delta = topValue
                self.selectedValues[rowID] = (topFamily, topValue, delta)
            else:
                self.selectedValues[rowID] = (None, None, None)

                # Method 1: Top Hit Above Threshold

    def selectTopAboveThreshold(self) -> None:
        """
        Select the top hit while using the family specific threshold cutoffs.
        """
        self.decoyDataframe = pd.read_csv(
            str(self.decoyStats),
            header=0,
            engine="c",
        )
        thresholdDict = dict(
            zip(self.decoyDataframe.Family, self.decoyDataframe[self.thresholdType])
        )

        self.selectedValues = {}
        filteredOutCount = 0

        for rowID, row in self.kmerCountTotals.iterrows():
            thresholdValues = row.index.map(thresholdDict.get)
            thresholdSeries = pd.Series(thresholdValues, index=row.index)

            rowValues = row[row > thresholdSeries]

            if not rowValues.empty:
                sortedRow = rowValues.sort_values(ascending=False)
                topValue = sortedRow.iloc[0]
                topFamily = sortedRow.index[0]
                if len(sortedRow) > 1:
                    secondValue = sortedRow.iloc[1]
                    delta = topValue - secondValue
                else:
                    delta = topValue - thresholdDict.get(topFamily, 0)
                self.selectedValues[rowID] = (topFamily, topValue, delta)
            else:
                self.selectedValues[rowID] = (None, None, None)
                filteredOutCount += 1

    # Method 2: Greatest Distance
    def selectByGreatestDistance(self) -> None:
        """
        Select the protein annotation by choosing the option that has the greatest distance from its family specific threshold.
        """
        self.decoyDataframe = pd.read_csv(
            str(self.decoyStats),
            header=0,
            engine="c",
        )
        thresholdDict = dict(
            zip(self.decoyDataframe.Family, self.decoyDataframe[self.thresholdType])
        )

        self.selectedValues = {}
        filteredOutCount = 0

        for rowID, row in self.kmerCountTotals.iterrows():
            distances = row - row.index.map(thresholdDict.get)
            positiveDistances = distances[distances > 0]

            if not positiveDistances.empty:
                greatestDistanceFamily = positiveDistances.idxmax()
                greatestDistanceValue = row[greatestDistanceFamily]
                delta = positiveDistances.max()
                self.selectedValues[rowID] = (
                    greatestDistanceFamily,
                    greatestDistanceValue,
                    delta,
                )
            else:
                self.selectedValues[rowID] = (None, None, None)
                filteredOutCount += 1

    
    # Method 3: Balanced Distance
    def selectByBalancedDistance(self, weightTop: float = 0.5, weightDistance: float = 0.5) -> None:
        """
        This method combines the previous two. In a weighted manner, choose the top protein annotation.
        Takes a combination of highest cosine similarity (top) and distance from family specific threshold.
        """
        self.decoyDataframe = pd.read_csv(
            str(self.decoyStats),
            header=0,
            engine="c",
        )
        thresholdDict = dict(
            zip(self.decoyDataframe.Family, self.decoyDataframe[self.thresholdType])
        )

        self.selectedValues = {}
        filteredOutCount = 0

        for rowID, row in self.kmerCountTotals.iterrows():
            thresholdValues = row.index.map(thresholdDict.get)
            thresholdSeries = pd.Series(thresholdValues, index=row.index)

            rowValuesAboveThreshold = row[row > thresholdSeries]
            distances = row - thresholdSeries
            positiveDistances = distances[distances > 0]

            candidates = {}

            if not rowValuesAboveThreshold.empty:
                topValue = rowValuesAboveThreshold.max()
                topFamily = rowValuesAboveThreshold.idxmax()
                top_threshold = thresholdDict.get(topFamily, 0)
                candidates[topFamily] = {
                    "value": topValue,
                    "delta": topValue - top_threshold,
                    "score": (topValue * weightTop)
                    + ((topValue - top_threshold) * weightDistance),
                }

            if not positiveDistances.empty:
                greatestDistanceFamily = positiveDistances.idxmax()
                greatestDistanceValue = row[greatestDistanceFamily]
                greatestDistanceThreshold = thresholdDict.get(
                    greatestDistanceFamily, 0
                )
                candidates[greatestDistanceFamily] = {
                    "value": greatestDistanceValue,
                    "delta": positiveDistances.max(),
                    "score": (greatestDistanceValue * weightTop)
                    + (positiveDistances.max() * weightDistance),
                }

            if candidates:
                bestCandidate = max(
                    candidates.items(), key=lambda x: x[1]["score"]
                )
                self.selectedValues[rowID] = (
                    bestCandidate[0],
                    bestCandidate[1]["value"],
                    bestCandidate[1]["delta"],
                )
            else:
                self.selectedValues[rowID] = (None, None, None)
                filteredOutCount += 1

    def formatAndWriteOutput(self) -> None:
        """
        Format and write the final output from Apply.
        """
        if config["learnapp"]["save_apply_associations"]:
            kmerCountTotals_write = pa.Table.from_pandas(
                self.kmerCountTotals
            )
            csv.write_csv(kmerCountTotals_write, self.outputSeqAnn)

        globalConfidenceScores = pd.read_csv(
            str(self.confidenceAssociations)
        )
        globalConfidenceScores.index = globalConfidenceScores[
            globalConfidenceScores.columns[0]
        ]
        globalConfidenceScores = globalConfidenceScores.iloc[:, 1:]
        globalConfidenceScores = globalConfidenceScores[
            globalConfidenceScores.columns[0]
        ].squeeze()

        resultsList = []
        for rowID in self.kmerCountTotals.index:
            if rowID in self.selectedValues:
                prediction, score, delta = self.selectedValues[rowID]
                if delta is None:
                    delta = 0
            else:
                prediction, score, delta = None, None, 0
            resultsList.append(
                {
                    "Sequence": rowID,
                    "Prediction": prediction,
                    "Score": score,
                    "delta": round(delta, 2),
                }
            )

        results = pd.DataFrame(resultsList)
        results.set_index("Sequence", inplace=True)

        results["Confidence"] = results["delta"].map(globalConfidenceScores)

        results.reset_index(inplace=True)
        resultsWrite = pa.Table.from_pandas(results)
        csv.write_csv(resultsWrite, self.outputKmerSummary)

    def execute_all(self) -> None:
        """
        Run all previously defined functions in sequence.
        """
        self.loadData()
        self.generateKmerCounts()
        self.constructKmerCountsDataframe()
        self.matchKmerCountsFormat()
        self.cosineSimilarity()
        # Select method based on selectionType
        if self.selectionType == "top_hit":
            if self.thresholdType is None:
                self.selectTopNoThreshold()
            else:
                self.selectTopAboveThreshold()
        elif self.selectionType == "greatest_distance":
            self.selectByGreatestDistance()
        elif self.selectionType == "combined_distance":
            weightTop = config["learnapp"].get("weight_top", 0.5)
            weightDistance = config["learnapp"].get("weight_distance", 0.5)
            self.selectByBalancedDistance(weightTop, weightDistance)
        else:
            raise ValueError(f"Invalid selectionType: {self.selectionType}")
        self.formatAndWriteOutput()


apply = KmerCompare(
    snakemake.input.compareAssociations,
    snakemake.input.data,
    snakemake.input.confidenceAssociations,
    snakemake.input.decoyStats,
    snakemake.input.annotation,
    snakemake.output.seqAnn,
    snakemake.output.kmerSummary,
    snakemake.params.selectionType,
    snakemake.params.thresholdType,
)
apply.execute_all()
