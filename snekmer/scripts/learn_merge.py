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


class Merge:
    """
Initializes the Merge object.

This object is designed to merge all dataframes containing kmer counts
and, if present, merge them with a base dataframe.

Attributes:
countsFiles (list): List of paths to the CSV files containing k-mer counts.
baseCountsPath (str): Path to the base CSV file for merging.
outputPath (str): Path to save the merged dataframe.
runningMerge (DataFrame or None): Merged dataframe from countsFiles.
baseCheck (bool): Flag to check if the base file exists and can be merged.
baseKmerCounts (DataFrame or None): Dataframe loaded from baseCountsPath.
"""

    def __init__(self, countsFiles, baseCountsPath, outputPath):
        self.countsFiles = countsFiles
        self.baseCountsPath = baseCountsPath
        self.outputPath = outputPath
        self.runningMerge = None
        self.baseCheck = False
        self.baseKmerCounts = None

    def mergeDataframes(self):
        """
Merges dataframes from the list of kmer count files.

Loads each dataframe from countsFiles, then successively merges them
into a running merged dataframe.
"""
        for fileNum, f in enumerate(self.countsFiles):
            kmerCounts = pd.read_csv(
                f,
                index_col="__index_level_0__",
                header=0,
                engine="pyarrow",
                na_values=[""],
            )
            kmerCounts.fillna(0, inplace=True)
            if fileNum == 0:
                self.runningMerge = kmerCounts
            else:
                self.runningMerge = (
                    pd.concat([self.runningMerge, kmerCounts])
                    .reset_index()
                    .groupby("__index_level_0__", sort=False)
                    .sum(min_count=1)
                ).fillna(0)
            print(
                f"Dataframes merged: {fileNum} out of {len(self.countsFiles)}"
            )

    def checkForBaseFile(self):
        """
Checks for the presence of a base file to merge with.

Sets the baseCheck flag to True if a CSV file is detected in the base path.
"""
        print("\nChecking for base file to merge with.\n")
        if "csv" in str(self.baseCountsPath):
            print(
                "CSV detected. Matching annotations, kmers, and totals will be summed. New annotations and kmers will be added."
            )
            self.baseCheck = True
        elif self.baseCountsPath == "":
            print("No base directory detected\n")
        elif str(self.baseCountsPath) == "input/base":
            print("Empty base directory detected\n")
        else:
            print(
                "No file type detected. Please use a .csv file in input/base directory.\n"
            )

    def confirmKmerCountsAndAlphabet(self):
        """
Confirms consistency between the alphabets and k-mer lengths
of the running merged dataframe and the base dataframe.

If any inconsistency is found, the baseCheck flag is set to False.
"""
        if self.baseCheck:
            self.baseKmerCounts = pd.read_csv(
                str(self.baseCountsPath),
                index_col="__index_level_0__",
                header=0,
                engine="pyarrow",
            )
            print("\nBase Database: \n")
            print(self.baseKmerCounts)
            check_1 = len(self.runningMerge.columns.values)
            alphabetInitial = set(
                itertools.chain(
                    *[
                        list(x)
                        for x in self.runningMerge.columns.values[3:check_1]
                    ]
                )
            )
            alphabet_base = set(
                itertools.chain(
                    *[
                        list(x)
                        for x in self.baseKmerCounts.columns.values[3:check_1]
                    ]
                )
            )
            if alphabet_base != alphabetInitial:
                self.baseCheck = False
                print("Different Alphabets Detected. Base File not merged.")
            if len(str(self.runningMerge.columns.values[1])) != len(
                str(self.baseKmerCounts.columns.values[1])
            ):
                self.baseCheck = False
                print("Different kmer lengths detected. Base File not merged.")

    def mergeWithBase(self):
        """
Merges the running merged dataframe with the base dataframe,
if the baseCheck flag is True.

If the flag is False, only the running merged dataframe is saved to output.
"""
        if self.baseCheck:
            print("\nMerged Database \n")
            xy = (
                pd.concat([self.baseKmerCounts, self.runningMerge])
                .reset_index()
                .groupby("__index_level_0__", sort=False)
                .sum(min_count=1)
            ).fillna(0)
            xy_out = pa.Table.from_pandas(xy, preserve_index=True)
            csv.write_csv(xy_out, self.outputPath)
        else:
            print("\nDatabase Merged. Not merged with base file.\n")
            runningMergeOut = pa.Table.from_pandas(
                self.runningMerge, preserve_index=True
            )
            csv.write_csv(runningMergeOut, self.outputPath)

    def executeAll(self):
        """
Executes all the merging steps in sequence.

This includes:
    1. Merging individual count dataframes.
    2. Checking for a base file.
    3. Confirming kmer counts and alphabet consistency.
    4. Merging with the base file if applicable.
"""
        self.mergeDataframes()
        self.checkForBaseFile()
        self.confirmKmerCountsAndAlphabet()
        self.mergeWithBase()


merger = Merge(snakemake.input.counts, snakemake.input.baseCounts, snakemake.output.totals)
merger.executeAll()