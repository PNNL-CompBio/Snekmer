
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
import re
import gzip
# ---------------------------------------------------------
# Files and Parameters
# ---------------------------------------------------------

config = snakemake.config

# change matplotlib backend to non-interactive
plt.switch_backend("Agg")

# ---------------------------------------------------------
# Run script
# ---------------------------------------------------------


class KmerCompare:
    """
Initializes the KmerCompare object.

This object is designed to compare kmer counts with provided annotations.

Attributes:
compareAssociations (str): Path to a CSV file containing kmer counts totals matrix.
annotationFiles (list): List of paths to files containing sequence annotations.
inputData (str): Path to input data for kmer analysis.
outputPath (str): Path to save the result.
annotation (list): List of dataframes loaded from annotationFiles.
kmerCountTotals (DataFrame or None): DataFrame of kmer count totals from compareAssociations.
seqAnnot (dict): Dictionary mapping sequence IDs to annotations.
kmerList (list): List of unique kmers found in input data.
seqKmerdict (dict): Dictionary mapping sequence IDs to their kmer counts.
totalSeqs (int): Total number of sequences processed.
kmerTotals (list): Total counts for each kmer across all sequences.
"""

    def __init__(
        self, compareAssociations, annotationFiles, inputData, outputPath
    ):
        self.compareAssociations = compareAssociations
        self.annotationFiles = annotationFiles
        self.inputData = inputData
        self.outputPath = outputPath
        self.annotation = []
        self.kmerCountTotals = None
        self.seqAnnot = {}
        self.kmerList = []
        self.seqKmerdict = {}
        self.totalSeqs = 0
        self.kmerTotals = []

    def generateInputs(self):
        """
Generates the necessary inputs for comparison.

Loads kmer counts and annotations into appropriate data structures.
"""
        self.kmerCountTotals = pd.read_csv(
            str(self.compareAssociations),
            index_col="__index_level_0__",
            header=0,
            engine="c",
        )
        for f in self.annotationFiles:
            self.annotation.append(pd.read_table(f))
        seqs = self.annotation[0]["id"].tolist()
        anns = self.annotation[0]["Family"].tolist()

        for i, seqid in enumerate(seqs):
            self.seqAnnot[seqid] = anns[i]

    def generateKmerCounts(self):
        """
Generates a dictionary of kmer counts for each sequence.

Processes the input data to count the occurrence of each kmer in each sequence.
"""
        kmerList, df = skm.io.load_npz(self.inputData)
        self.kmerList = kmerList[0]
        seqids = df["sequence_id"]
        self.kmerTotals = [0] * len(self.kmerList)
        k_len = len(self.kmerList[0])

        for i, seq in enumerate(seqids):
            v = df["sequence"][i]
            kCounts = {}
            items = [
                v[item : (item + k_len)]
                for item in range(0, (len((v)) - k_len + 1))
            ]
            for j in items:
                kCounts[j] = kCounts.get(j, 0) + 1
            store = [
                kCounts[item] if item in kCounts else 0
                for item in self.kmerList
            ]
            for i, item in enumerate(self.kmerList):
                if item in kCounts:
                    self.kmerTotals[i] += kCounts[item]
            self.seqKmerdict[seq] = store

    def addKnownTag(self):
        """
Modifies sequence IDs to include a known or unknown tag based on annotation.

Uses annotations to determine whether a sequence is known or unknown.
Appends a tag to the sequence ID accordingly.
"""
        seqs = set(self.annotation[0]["id"].tolist())
        count = 0
        for seqid in list(self.seqKmerdict):
            x = re.findall(r"\|(.*?)\|", seqid)[0]
            if x not in seqs:
                self.seqKmerdict[(str(x) + "_unknown_" + str(count))] = (
                    self.seqKmerdict.pop(seqid)
                )
            else:
                self.seqKmerdict[
                    (str(self.seqAnnot[x]) + "_known_" + str(count))
                ] = self.seqKmerdict.pop(seqid)
            count += 1
        self.totalSeqs = len(self.seqKmerdict)

    def constructKmerCountsDataframe(self):
        """
Constructs a pandas DataFrame of kmer counts for each sequence.

Returns:
    DataFrame: A DataFrame where rows represent sequences (and a total row),
            and columns represent kmers.
"""
        kmerCounts = pd.DataFrame(self.seqKmerdict.values())
        kmerCounts.insert(0, "Annotations", 1, True)
        self.kmerTotals.insert(0, self.totalSeqs)
        kmerCounts = pd.DataFrame(
            np.insert(kmerCounts.values, 0, values=self.kmerTotals, axis=0)
        )
        kmerCounts.columns = ["Sequence count"] + list(self.kmerList)
        kmerCounts.index = ["Totals"] + list(self.seqKmerdict.keys())
        return kmerCounts

    def matchkmerCountsFormat(self, kmerCounts):
        """
Matches the format of the provided kmer counts DataFrame to the format of the
comparison data. Ensures columns align correctly.

Args:
    kmerCounts (DataFrame): DataFrame of kmer counts to format.

Returns:
    DataFrame: Formatted kmer counts DataFrame.
"""
        if len(str(kmerCounts.columns.values[10])) == len(
            str(self.kmerCountTotals.columns.values[10])
        ):
            compareCheck = True
        else:
            compareCheck = False

        if compareCheck:
            check_1 = len(kmerCounts.columns.values)
            alphabetInitial = set(
                itertools.chain(
                    *[list(x) for x in kmerCounts.columns.values[10:check_1]]
                )
            )
            alphabetCompare = set(
                itertools.chain(
                    *[
                        list(x)
                        for x in self.kmerCountTotals.columns.values[10:check_1]
                    ]
                )
            )
            if alphabetCompare != alphabetInitial:
                compareCheck = False

        if not compareCheck:
            print("Compare Check Failed. ")
            sys.exit()

        kmerCounts.drop("Totals", axis=0, inplace=True)
        kmerCounts.drop("Sequence count", axis=1, inplace=True)

        self.kmerCountTotals.drop("Totals", axis=0, inplace=True)
        self.kmerCountTotals.drop("Kmer Count", axis=1, inplace=True)
        self.kmerCountTotals.drop("Sequence count", axis=1, inplace=True)

        columnOrder = list(
            set(kmerCounts.columns) | set(self.kmerCountTotals.columns)
        )
        kmerCounts = kmerCounts.reindex(columns=columnOrder, fill_value=0)
        self.kmerCountTotals = self.kmerCountTotals.reindex(
            columns=columnOrder, fill_value=0
        )

        return kmerCounts

    def calculateCosineSimilarity(self, kmerCounts):
        """
Calculates the cosine similarity between the input kmer counts and comparison data.

Args:
    kmerCounts (DataFrame): DataFrame of kmer counts for comparison.

Returns:
    DataFrame: A DataFrame of cosine similarity scores.
"""
        cosine_df = sklearn.metrics.pairwise.cosine_similarity(
            self.kmerCountTotals, kmerCounts
        ).T
        finalMatrixWithScores = pd.DataFrame(
            cosine_df,
            columns=self.kmerCountTotals.index,
            index=kmerCounts.index,
        )
        return finalMatrixWithScores

    def filterTopTwoValues(self, finalMatrixWithScores):
        """
Filters the similarity scores to keep only the top two values for each row.

Args:
    finalMatrixWithScores (DataFrame): DataFrame of similarity scores.

Returns:
    DataFrame: DataFrame with all but the top two scores set to NaN.
"""
        topTwoIndices = np.argsort(-finalMatrixWithScores.values, axis=1)[:, :2]
        mask = np.zeros_like(finalMatrixWithScores.values, dtype=bool)
        for i, (indexOne, indexTwo) in enumerate(topTwoIndices):
            mask[i, indexOne] = True
            mask[i, indexTwo] = True
        finalMatrixWithScores.values[~mask] = np.nan
        return finalMatrixWithScores

    def filter_top_one_value(self, finalMatrixWithScores):
        """
Filters the similarity scores to keep only the top one value for each row.

Args:
    finalMatrixWithScores (DataFrame): DataFrame of similarity scores.

Returns:
    DataFrame: DataFrame with all but the top single score set to NaN.
"""
        top_1_indices = np.argsort(-finalMatrixWithScores.values, axis=1)[:, :1]
        mask = np.zeros_like(finalMatrixWithScores.values, dtype=bool)

        for i, indexOne in enumerate(top_1_indices):
            mask[i, indexOne] = True

        finalMatrixWithScores.values[~mask] = np.nan
        return finalMatrixWithScores

    def writeOutput(self, finalMatrixWithScores):
        """
Writes the provided DataFrame to a CSV file at the specified output path.

Args:
    finalMatrixWithScores (DataFrame): DataFrame to write to CSV.
"""
        finalMatrixWithScoresWrite = pa.Table.from_pandas(finalMatrixWithScores)
        with gzip.open(self.outputPath, "wb") as gzipped_file:
            csv.write_csv(finalMatrixWithScoresWrite, gzipped_file)

    def executeAll(self, config):
        """
Executes all the comparison steps in sequence.

This includes:
    1. Loading inputs.
    2. Generating kmer counts.
    # 3. Adding known/unknown tags.
    4. Constructing a k-mer counts dataframe.
    5. Matching format with comparison data.
    6. Calculating cosine similarity.
    7. Filtering to keep top two values (if applicable).
    8. Writing results to output.
"""
        self.generateInputs()
        self.generateKmerCounts()
        self.addKnownTag()
        kmerCounts = self.constructKmerCountsDataframe()
        kmerCounts = self.matchkmerCountsFormat(kmerCounts)
        finalMatrixWithScores = self.calculateCosineSimilarity(kmerCounts)
        if not config["learnapp"]["save_apply_associations"]:
            finalMatrixWithScores = self.filterTopTwoValues(
                finalMatrixWithScores
            )
        self.writeOutput(finalMatrixWithScores)


analysis = KmerCompare(
    snakemake.input.compareAssociations, snakemake.input.annotation, snakemake.input.data, snakemake.output.apply
)
analysis.executeAll(config)
