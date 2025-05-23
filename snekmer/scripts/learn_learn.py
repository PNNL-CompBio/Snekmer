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
# ---------------------------------------------------------
# Files and Parameters
# ---------------------------------------------------------

config = snakemake.config

# change matplotlib backend to non-interactive
plt.switch_backend("Agg")

outDir = skm.io.define_output_dir(
    config["alphabet"], config["k"], nested=config["nested_output"]
)

# ---------------------------------------------------------
# Run script
# ---------------------------------------------------------

class Library:
    """
Initializes a Library object.

This object is responsible for processing annotations and kmer data.
It creates a kmer-count matrix for each input file.

Attributes:
annotation (list): A list to store annotations loaded from files.
seqAnnot (dict): A dictionary mapping sequence IDs to their annotations.
kmerList (list): A list of unique kmers present in the data.
df (DataFrame or None): DataFrame containing kmer data.
seqids (list): A list of sequence IDs.
kmerTotals (list): A list to store total counts of each k-mer across all sequences.
seqKmerdict (dict): A dictionary mapping sequence IDs to their k-mer counts.
annotationCounts (dict): A dictionary mapping annotations to their counts.
totalSeqs (int): The total number of sequences after filtering.
"""

    def __init__(self):
        self.annotation = []
        self.seqAnnot = {}
        self.kmerList = []
        self.df = None
        self.seqids = []
        self.kmerTotals = []
        self.seqKmerdict = {}
        self.annotationCounts = {}
        self.totalSeqs = 0

    def loadAnnotations(self, input_annotation):
        """
Load annotations from a list of provided input files.

Args:
    input_annotation (list): List of file paths containing annotations.
"""
        for f in input_annotation:
            self.annotation.append(pd.read_table(f))
        annotations = pd.concat(self.annotation)
        seqs = annotations["id"].tolist()
        # anns = annotations["Family"].tolist()
        anns = [str(fam) for fam in annotations["Family"].tolist()]

        for i, seqid in enumerate(seqs):
            self.seqAnnot[seqid] = anns[i]
        self.seqs = set(seqs)

    def loadData(self, inputData):
        """
Load and format kmer data from the provided input.

Args:
    inputData (str): Path to the data file.
"""
        self.kmerList, self.df = skm.io.load_npz(inputData)
        self.kmerList = self.kmerList[0]
        self.seqids = self.df["sequence_id"]
        for item in self.kmerList:
            self.kmerTotals.append(0)

    def generateKmerCounts(self):
        """
Generate kmer counts for sequences present in the data.
"""
        k_len = len(self.kmerList[0])
        for i, seq in enumerate(self.seqids):
            v = self.df["sequence"][i]
            kCounts = self._computeKmerCountsForSequence(v, k_len)
            self.seqKmerdict[seq] = kCounts

    def filterAndConstruct(self):
        """
Filters sequences not present in annotations and constructs annotation counts.
"""
        self.totalSeqs = len(self.seqKmerdict)
        for i, seqid in enumerate(list(self.seqKmerdict)):
            x = re.findall(r"\|(.*?)\|", seqid)[
                0
            ]  # A Note, it could be useful to allow a user to define their own regex.
            if x not in self.seqs:
                del self.seqKmerdict[seqid]
            else:
                self._processAnnotationCounts(seqid, x)

    def formatAndWriteOutput(self, inputData):
        """
Writes processed kmer counts to an output CSV file.

Args:
    inputData (str): Path to the data file (used for naming the output file).
"""
        kmerCounts = pd.DataFrame(self.seqKmerdict.values())
        kmerCounts.insert(
            0, "Annotations", self.annotationCounts.values(), True
        )

        kmerCountsValues = (
            kmerCounts[list(kmerCounts.columns[1:])].sum(axis=1).to_list()
        )
        kmerCounts.insert(1, "Kmer Count", kmerCountsValues, True)

        self.kmerTotals[0:0] = [self.totalSeqs, sum(self.kmerTotals)]
        colnames = ["Sequence count"] + ["Kmer Count"] + list(self.kmerList)
        kmerCounts = pd.DataFrame(
            np.insert(kmerCounts.values, 0, values=self.kmerTotals, axis=0)
        )
        kmerCounts.columns = colnames
        newIndex = ["Totals"] + list(self.annotationCounts.keys())
        kmerCounts.index = newIndex
        kmerCounts.replace(0, "", inplace=True)
        out_name = join(
            outDir, "learn", "kmer-counts-" + str(inputData)[14:-4] + ".csv"
        )
        kmerCounts.index.name = "__index_level_0__"
        kmerCounts.to_csv(out_name, index=True)

    def _computeKmerCountsForSequence(self, v, k_len):
        """
Computes k-mer counts for a given sequence.

Args:
    v (str): The sequence.
    k_len (int): Length of the k-mer.

Returns:
    list: List of k-mer counts for the sequence.
"""
        items = [
            v[item : item + k_len] for item in range(0, len(v) - k_len + 1)
        ]
        kCounts = {}
        for j in items:
            kCounts[j] = kCounts.get(j, 0) + 1
        store = []
        for i, item in enumerate(self.kmerList):
            if item in kCounts:
                store.append(kCounts[item])
                self.kmerTotals[i] += kCounts[item]
            else:
                store.append(0)
        return store

    def _processAnnotationCounts(self, seqid, x):
        """
Processes annotation counts by aggregating them based on annotation labels.

Args:
    seqid (str): Sequence ID.
    x (str): Extracted annotation ID from seqid.
"""
        if self.seqAnnot[x] not in self.seqKmerdict:
            self.seqKmerdict[self.seqAnnot[x]] = self.seqKmerdict.pop(seqid)
        else:
            zipped_lists = zip(
                self.seqKmerdict.pop(seqid),
                self.seqKmerdict[self.seqAnnot[x]],
            )
            self.seqKmerdict[self.seqAnnot[x]] = [
                sum(pair) for pair in zipped_lists
            ]
        if self.seqAnnot[x] not in self.annotationCounts:
            self.annotationCounts[self.seqAnnot[x]] = 1
        else:
            self.annotationCounts[self.seqAnnot[x]] += 1

    def executeAll(self, input_annotation, inputData):
        """
Execute the entire sequence of operations in the Library process.

Args:
    input_annotation (list): List of file paths containing annotations.
    inputData (str): Path to the data file.
"""
        self.loadAnnotations(input_annotation)
        self.loadData(inputData)
        self.generateKmerCounts()
        self.filterAndConstruct()
        self.formatAndWriteOutput(inputData)


library = Library()
library.executeAll(snakemake.input.annotation, snakemake.input.data)
