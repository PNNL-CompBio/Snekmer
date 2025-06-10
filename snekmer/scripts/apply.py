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

from snekmer.learn import generateKmerCounts, matchKmerCountsFormat, applySelectionMethod

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
        self.kmerTotals = None
        self.seqKmerDict = {}

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


    def build_threshold_dict(self):
        if self.thresholdType is None:
            self.thresholdDict = {}
        else:
            df = pd.read_csv(str(self.decoyStats), engine="c")
            self.thresholdDict = dict(zip(df.Family.astype(str),
                                        df[self.thresholdType]))
            
                        
    def execute_all(self) -> None:
        weightTop = config["learnapp"].get("weight_top", 0.5)
        weightDistance = config["learnapp"].get("weight_distance", 0.5)
        self.loadData()
        
        self.kmerList, self.kmerTotals, self.seqKmerDict = generateKmerCounts(self.data, 
                                                                            self.kmerList,
                                                                            self.kmerTotals, 
                                                                            self.seqKmerDict, 
                                                                            False)
        self.constructKmerCountsDataframe()
        self.kmerCounts, self.kmerCountTotals = matchKmerCountsFormat(
            self.kmerCounts, self.kmerCountTotals)
        self.cosineSimilarity()
        self.build_threshold_dict()
        
        preds, deltas, topTwo = applySelectionMethod(
            self.kmerCountTotals,
            self.selectionType,
            self.thresholdType,
            self.thresholdDict,
            weightTop,
            weightDistance)
        
        
        topTwo.index = self.kmerCountTotals.index
        self.selectedValues = {}
        for seq_id, pred, delta in zip(self.kmerCountTotals.index, preds, deltas):
            if pred is None:
                self.selectedValues[seq_id] = (None, None, None)
            else:
                top_val = float(topTwo.loc[seq_id, "keyValueOne"])
                self.selectedValues[seq_id] = (pred, top_val, float(delta))

        self.formatAndWriteOutput()


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
