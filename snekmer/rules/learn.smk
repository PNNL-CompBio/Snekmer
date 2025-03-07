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
from itertools import product, repeat
from multiprocessing import Pool
from os import makedirs
from os.path import basename, dirname, exists, join, split, splitext
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

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


# Note:
# Pyarrow installed via "conda install -c conda-forge pyarrow"
# collect all fasta-like files, unzipped filenames, and basenames
inputDir = (
    "input"
    if (("input_dir" not in config) or (str(config["input_dir"]) == "None"))
    else config["input_dir"]
)
inputFiles = glob(join(inputDir, "*"))
# base_file = glob(join(inputDir,"base" "*"))
zipped = [fa for fa in inputFiles if fa.endswith(".gz")]
unzipped = [
    fa.rstrip(".gz")
    for fa, ext in product(inputFiles, config["input_file_exts"])
    if fa.rstrip(".gz").endswith(f".{ext}")
]
annotFiles = glob(join("annotations", "*.ann"))
baseCounts = glob(join("base", "counts", "*.csv"))
baseConfidence = glob(join("base", "confidence", "*.csv"))
baseFamilyCheckpoint = glob(join("base", "thresholds", "*.csv"))

uzMap = {
    skm.utils.split_file_ext(f)[0]: skm.utils.split_file_ext(f)[1] for f in zipped
}
faMap = {
    skm.utils.split_file_ext(f)[0]: skm.utils.split_file_ext(f)[1] for f in unzipped
}

# get unzipped filenames
UZS = [f"{f}.{ext}" for f, ext in uzMap.items()]
# isolate basenames for all files
FAS = list(faMap.keys())
# parse any background files
backgroundFiles = glob(join(inputDir, "background", "*"))
if len(backgroundFiles) > 0:
    backgroundFiles = [skm.utils.split_file_ext(basename(f))[0] for f in backgroundFiles]
# terminate with error if invalid alphabet specified
skm.alphabet.check_valid(config["alphabet"])
# define output directory (helpful for multiple runs)
outDir = skm.io.define_output_dir(
    config["alphabet"], config["k"], nested=config["nested_output"]
)

output_prefixes = (
    ["vector"] if not config["learnapp"]["fragmentation"] else ["vector", "vector_frag"]
)

if config["learnapp"]["selection"] != "top_hit" and config["learnapp"]["threshold"] == "None":
    raise Exception("The only selection method that allows for None is `top_hit`. Other methods inherently use a threshold")

# define output files to be created by snekmer
rule all:
    input:
        # Vector files
        expand(join(outDir, "vector", "{nb}.npz"), nb=FAS),
        # Kmer counts files for learning
        expand(join(outDir, "learn", "kmer-counts-{nb}.csv"), nb=FAS),
        join(outDir, "learn", "kmer-counts-total.csv"),
        # Fragmentation outputs (only if enabled)
        expand(join(outDir, "fragmented", "{nb}.fasta"), nb=FAS)
            if config["learnapp"]["fragmentation"] else [],
        expand(join(outDir, "vector_frag", "{nb}.npz"), nb=FAS)
            if config["learnapp"]["fragmentation"] else [],
        # Evaluation application sequences (or fragments if fragmentation is enabled)
        expand(
            join(
                outDir,
                "eval_apply_sequences" if not config["learnapp"]["fragmentation"] else "eval_apply_frag",
                "seq-annotation-scores-{nb}.csv.gz"
            ),
            nb=FAS,
        ),
        # Evaluation configuration outputs
        join("output", "eval_conf", "family_summary_stats.csv"),
        join("output", "eval_conf", "global-confidence-scores.csv"),
        # Duplicated some outputs for each of use on the user side.
        join(outDir, "apply_inputs", "counts", "kmer-counts-total.csv"),
        join(outDir, "apply_inputs", "stats", "family_summary_stats.csv"),
        join(outDir, "apply_inputs", "confidence", "global-confidence-scores.csv")


# if any files are gzip zipped, unzip them
use rule unzip from process with:
    output:
        unzipped=join(inputDir, "{uz}"),
        zipped=join(inputDir, "zipped", "{uz}.gz"),


if config["learnapp"]["fragmentation"]:

    rule fragmentation:
        input:
            fasta=lambda wildcards: join(
                inputDir, f"{wildcards.nb}.{faMap[wildcards.nb]}"
            ),
        output:
            fasta_out=join(outDir, "fragmented", "{nb}.fasta"),
        params:
            version=config["learnapp"]["version"],
            frag_length=config["learnapp"]["frag_length"],
            location=config["learnapp"]["location"],
            min_length=config["learnapp"]["min_length"],
            seed=config["learnapp"]["seed"],
        run:
            random.seed(params.seed)  # setting the seed for randomness

            def fragment(sequence, version, length, location, min_length):
                """
                Fragment a given sequence.

                Parameters:
                - sequence (str): the sequence to fragment.
                - version (str): fragmentation method, either "absolute" or "percent".
                - length (int): length for fragmentation. If version is "percent", this is treated as a percentage.
                - location (str): where the fragmentation happens - "start", "end", or "random".
                - min_length (int): minimum length a fragment should be to be retained.

                Returns:
                - list of fragments.
                """
                frags = []
                if version == "absolute":
                    for i in range(0, len(sequence) - length + 1):
                        frags.append(sequence[i : i + length])

                elif version == "percent":
                    actual_length = int(len(sequence) * (length / 100))
                    for i in range(0, len(sequence) - actual_length + 1):
                        frags.append(sequence[i : i + actual_length])

                        # Filter fragments based on min_length
                frags = [frag for frag in frags if len(frag) >= min_length]

                # Retention logic based on the location parameter
                if location == "start":
                    if frags:
                        frags = [frags[0]]
                elif location == "end":
                    if frags:
                        frags = [frags[-1]]
                elif location == "random":
                    if frags:
                        chosenIndex = random.randint(0, len(frags) - 1)
                        frags = [frags[chosenIndex]]

                return frags

            with open(input.fasta, "r") as f:
                fastaSequences = SeqIO.parse(f, "fasta")

                with open(output.fasta_out, "w") as the_file:
                    for fasta in fastaSequences:
                        title_line, sequence = fasta.description, str(fasta.seq)

                        fragments = fragment(
                            sequence,
                            params.version,
                            params.frag_length,
                            params.location,
                            params.min_length,
                        )

                        for i, frag in enumerate(fragments):
                            the_file.write(
                                ">" + title_line + " Fragment=" + str(i) + "\n"
                            )
                            the_file.write(frag + "\n")


use rule vectorize from kmerize with:
    input:
        fasta=lambda wildcards: join(
            outDir if wildcards.prefix == "vector_frag" else inputDir,
            "fragmented" if wildcards.prefix == "vector_frag" else "",
            f"{wildcards.nb}.{faMap[wildcards.nb]}",
        ),
    output:
        data=join(outDir, "{prefix}", "{nb}.npz"),
        kmerobj=join(outDir, "kmerize_{prefix}", "{nb}.kmers"),
    log:
        join(outDir, "{prefix}_kmerize", "log", "{nb}.log"),

# WORKFLOW to learn kmer associations
rule learn:
    input:
        data = join(outDir, "vector", "{nb}.npz"),
        annotation=expand("{an}", an=annotFiles),
    output:
        counts = join(outDir, "learn", "kmer-counts-{nb}.csv"),
    log:
        join(outDir, "learn", "log", "learn-{nb}.log"),
    run:
        start_time = datetime.now()
        with open(log[0], "a") as f:
            f.write(f"start time:\t{start_time}\n")


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
                    x = re.findall(r"\|(.*?)\|", seqid)[0]      # A Note, it could be useful to allow a user to define their own regex.
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
                out_name = join(outDir, "learn", "kmer-counts-" + str(inputData)[14:-4] + ".csv")
                kmerCounts.index.name = "__index_level_0__"
                kmerCounts.to_csv(out_name, index=True)
                skm.utils.log_runtime(log[0], start_time)

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
                    self.seqKmerdict[self.seqAnnot[x]] = self.seqKmerdict.pop(
                        seqid
                    )
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
        library.executeAll(input.annotation, input.data)


rule merge:
    input:
        counts=expand(join(outDir, "learn", "kmer-counts-{nb}.csv"), nb=FAS),
        baseCounts=expand("{bf}", bf=baseCounts),
    output:
        totals=join(outDir, "learn", "kmer-counts-total.csv"),
    log:
        join(outDir, "learn", "log", "merge.log"),
    run:
        start_time = datetime.now()
        with open(log[0], "a") as f:
            f.write(f"start time:\t{start_time}\n")


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
                            *[list(x) for x in self.baseKmerCounts.columns.values[3:check_1]]
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


        merger = Merge(input.counts, input.baseCounts, output.totals)
        merger.executeAll()


rule eval_apply_reverse_seqs:
    input:
        data=join(
            outDir,
            "vector"
            if config["learnapp"]["fragmentation"] == False
            else "vector_frag",
            "{nb}.npz",
        ),
        annotation=expand("{an}", an=annotFiles),
        compareAssociations=join(outDir, "learn", "kmer-counts-total.csv"),

    output:
        apply=join(
            outDir,
            "eval_apply_reversed"
            if config["learnapp"]["fragmentation"] == False
            else "eval_apply_frag",
            "seq-annotation-scores-{nb}.csv.gz",
        ),
    log:
        join(outDir, "eval_apply_reversed", "log", "{nb}.log"),
    run:
        start_time = datetime.now()
        with open(log[0], "a") as f:
            f.write(f"start time:\t{start_time}\n")


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
                    

            def generate_reverse_kmerCounts(self):
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
                    v = v[::-1]
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
                                for x in self.kmerCountTotals.columns.values[
                                    10:check_1
                                ]
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
                finalMatrixWithScores = finalMatrixWithScores.round(3)  #Adding rounding for storage space.
                return finalMatrixWithScores

            def writeOutput(self, finalMatrixWithScores):
                """
        Writes the provided DataFrame to a CSV file at the specified output path.

        Args:
            finalMatrixWithScores (DataFrame): DataFrame to write to CSV.
        """
                finalMatrixWithScoresWrite = pa.Table.from_pandas(
                    finalMatrixWithScores
                )
                with gzip.open(self.outputPath, 'wb') as gzipped_file:
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
                self.generate_reverse_kmerCounts()
                kmerCounts = self.constructKmerCountsDataframe()
                kmerCounts = self.matchkmerCountsFormat(kmerCounts)
                finalMatrixWithScores = self.calculateCosineSimilarity(kmerCounts)
                self.writeOutput(finalMatrixWithScores)
                
        analysis = KmerCompare(
            input.compareAssociations, input.annotation, input.data, output.apply
        )
        analysis.executeAll(config)
        skm.utils.log_runtime(log[0], start_time)

rule reverseDecoy_evaluations:
    input:
        evalApplyData=expand(
            join(
                outDir,
                "eval_apply_reversed"
                if config["learnapp"]["fragmentation"] == False
                else "eval_apply_frag",
                "seq-annotation-scores-{nb}.csv.gz",
            ),
            nb=FAS,
        ),
        baseFamilyCheckpoint=expand("{bc}", bc=baseFamilyCheckpoint),
    output:
        familyStats = join(outDir, "eval_conf", "family_summary_stats.csv"),
        checkpoint = join(outDir, "eval_conf", "family_stats_checkpoint.csv")
    run:

        def load_existing_stats_from_csv(csv_file):
            """
            Load existing family statistics from a CSV checkpoint file.

            The CSV must have the following columns:
            Family, count, sum, sumSqr, min, max, values_for_percentiles (JSON-encoded list)
            """
            csv_file = str(csv_file)
            print(f"csv_file:{csv_file}")
            print(f"csv_file type:{type(csv_file)}")
            if not csv_file :
                return {}

            df = pd.read_csv(csv_file)
            combinedStats = {}
            for _, row in df.iterrows():
                family = row["Family"]
                combinedStats[family] = {
                    "count": int(row["count"]),
                    "sum": float(row["sum"]),
                    "sumSqr": float(row["sumSqr"]),
                    "min": float(row["min"]),
                    "max": float(row["max"]),
                    "values_for_percentiles": json.loads(row["values_for_percentiles"])
                }
            return combinedStats

        def save_stats_to_csv(combinedStats, csv_file):
            """
            Save the combined statistics to a CSV file.

            Each row will contain:
            Family, count, sum, sumSqr, min, max, values_for_percentiles (as a JSON string)
            """
            rows = []
            for family, stats in combinedStats.items():
                row = {
                    "Family": family,
                    "count": stats["count"],
                    "sum": stats["sum"],
                    "sumSqr": stats["sumSqr"],
                    "min": stats["min"],
                    "max": stats["max"],
                    # Encode values_for_percentiles as JSON string
                    "values_for_percentiles": json.dumps(stats["values_for_percentiles"])
                }
                rows.append(row)

            df = pd.DataFrame(rows)
            df.to_csv(csv_file, index=False)

        def collectFamilyStatistics(filename, existingStats=None):
            """
            Reads a CSV file in chunks and updates statistics for each family (column).
            """
            chunk_size = 10000  # Adjust based on available memory
            reservoir_size = 100000  # Size of the reservoir for percentiles

            if existingStats is None:
                existingStats = {}

            for chunk in pd.read_csv(filename, chunksize=chunk_size, engine="c"):
                families = chunk.columns[:-1]  # Exclude the last column if it's the sequence name

                for family in families:
                    values = chunk[family].dropna().astype(float).values
                    if family not in existingStats:
                        existingStats[family] = {
                            'count': 0,
                            'sum': 0.0,
                            'sumSqr': 0.0,
                            'min': np.inf,
                            'max': -np.inf,
                            'values_for_percentiles': []
                        }

                    stats = existingStats[family]

                    n = len(values)
                    if n == 0:
                        continue

                    # Update count and sum statistics
                    stats['sum'] += values.sum()
                    stats['sumSqr'] += np.dot(values, values)
                    stats['min'] = min(stats['min'], values.min())
                    stats['max'] = max(stats['max'], values.max())  

                    # Reservoir sampling for percentiles
                    for value in values:
                        stats['count'] += 1  # Update total count
                        total_seen = stats['count']

                        if len(stats['values_for_percentiles']) < reservoir_size:
                            # Fill the reservoir until it reaches the desired size
                            stats['values_for_percentiles'].append(value)
                        else:
                            # Replace elements with decreasing probability
                            j = random.randint(0, total_seen - 1)
                            if j < reservoir_size:
                                stats['values_for_percentiles'][j] = value

                del chunk  # Free memory

            return existingStats

        def generateFamilyStatistics(combinedStats):
            """
            Generates statistics for each family using the combined statistics.
            """
            statsData = {
                'Family': [],
                'Mean': [],
                'Std Dev': [],
                'Min': [],
                '10th Percentile': [],
                '20th Percentile': [],
                '25th Percentile': [],
                '30th Percentile': [],
                '40th Percentile': [],
                'Median': [],
                '60th Percentile': [],
                '70th Percentile': [],
                '75th Percentile': [],
                '80th Percentile': [],
                '90th Percentile': [],
                'Max': [],
                '1 Std Dev Above': [],
                '1 Std Dev Below': [],
                '2 Std Dev Above': [],
                '2 Std Dev Below': [],
            }

            for family, stats in combinedStats.items():
                n = stats['count']
                sum_ = stats['sum']
                sumSqr = stats['sumSqr']
                mean = sum_ / n if n > 0 else 0.0
                variance = (sumSqr - (sum_ ** 2) / n) / (n - 1) if n > 1 else 0.0
                std_dev = np.sqrt(variance)

                values = np.array(stats['values_for_percentiles'])
                if len(values) > 0:
                    percentiles = np.percentile(values, [10, 20, 25, 30, 40, 50, 60, 70, 75, 80, 90])
                else:
                    # If no values, fill with NaN
                    percentiles = [np.nan] * 11

                statsData['Family'].append(family)
                statsData['Mean'].append(round(mean, 3))
                statsData['Std Dev'].append(round(std_dev, 3))
                statsData['Min'].append(round(stats['min'], 3))
                statsData['10th Percentile'].append(round(percentiles[0], 3))
                statsData['20th Percentile'].append(round(percentiles[1], 3))
                statsData['25th Percentile'].append(round(percentiles[2], 3))
                statsData['30th Percentile'].append(round(percentiles[3], 3))
                statsData['40th Percentile'].append(round(percentiles[4], 3))
                statsData['Median'].append(round(percentiles[5], 3))
                statsData['60th Percentile'].append(round(percentiles[6], 3))
                statsData['70th Percentile'].append(round(percentiles[7], 3))
                statsData['75th Percentile'].append(round(percentiles[8], 3))
                statsData['80th Percentile'].append(round(percentiles[9], 3))
                statsData['90th Percentile'].append(round(percentiles[10], 3))
                statsData['Max'].append(round(stats['max'], 3))
                statsData['1 Std Dev Above'].append(round(mean + std_dev, 3))
                statsData['1 Std Dev Below'].append(round(mean - std_dev, 3))
                statsData['2 Std Dev Above'].append(round(mean + 2 * std_dev, 3))
                statsData['2 Std Dev Below'].append(round(mean - 2 * std_dev, 3))

            return pd.DataFrame(statsData)

        # Load existing stats from checkpoint CSV if exists
        if input.baseFamilyCheckpoint:
            print(f"input.baseFamilyCheckpoint is: {input.baseFamilyCheckpoint}")
            combinedStats = load_existing_stats_from_csv(input.baseFamilyCheckpoint)
        else:
            combinedStats = None

        # Update combinedStats with new data
        for filename in input.evalApplyData:
            combinedStats = collectFamilyStatistics(filename, existingStats=combinedStats)

        # Generate updated family statistics
        familyStatisticsDf = generateFamilyStatistics(combinedStats)
        familyStatisticsDf.to_csv(output.familyStats, index=False)

        # Save updated combinedStats to checkpoint CSV
        save_stats_to_csv(combinedStats, output.checkpoint)


rule eval_apply_sequences:
    input:
        data=join(
            outDir,
            "vector"
            if config["learnapp"]["fragmentation"] == False
            else "vector_frag",
            "{nb}.npz",
        ),
        annotation=expand("{an}", an=annotFiles),
        compareAssociations=join(outDir, "learn", "kmer-counts-total.csv"),
    output:
        apply=join(
            outDir,
            "eval_apply_sequences"
            if config["learnapp"]["fragmentation"] == False
            else "eval_apply_frag",
            "seq-annotation-scores-{nb}.csv.gz",
        ),
    log:
        join(outDir, "eval_apply_sequences", "log", "{nb}.log"),
    run:
        start_time = datetime.now()
        with open(log[0], "a") as f:
            f.write(f"start time:\t{start_time}\n")


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
                        self.seqKmerdict[
                            (str(x) + "_unknown_" + str(count))
                        ] = self.seqKmerdict.pop(seqid)
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
                                for x in self.kmerCountTotals.columns.values[
                                    10:check_1
                                ]
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
                topTwoIndices = np.argsort(-finalMatrixWithScores.values, axis=1)[
                    :, :2
                ]
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
                finalMatrixWithScoresWrite = pa.Table.from_pandas(
                    finalMatrixWithScores
                )
                with gzip.open(self.outputPath, 'wb') as gzipped_file:
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
            input.compareAssociations, input.annotation, input.data, output.apply
        )
        analysis.executeAll(config)
        skm.utils.log_runtime(log[0], start_time)


rule evaluate:
    input:
        evalApplyData=expand(
            join(
                outDir,
                "eval_apply_sequences"
                if config["learnapp"]["fragmentation"] == False
                else "eval_apply_frag",
                "seq-annotation-scores-{nb}.csv.gz",
            ),
            nb=FAS,
        ),
        baseConfidence=expand("{bc}", bc=baseConfidence),
        reverseDecoyStats = join(outDir, "eval_conf", "family_summary_stats.csv"),
    output:
        evalGlob = join(outDir, "eval_conf", "global-confidence-scores.csv"),
    params:
        modifier = config["learnapp"]["conf_weight_modifier"],
    log:
        join(outDir, "eval_conf", "log", "conf.log"),
    run:
        start_time = datetime.now()
        with open(log[0], "a") as f:
            f.write(f"start time:\t{start_time}\n")

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
                """
                Initializes the Evaluator object with input data and paths.
                """
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
                    threshold_dict = dict(zip(thresholds_df.Family.astype(str), thresholds_df[threshold_type]))

                    self.threshold_dict = threshold_dict  # Store for later use

                # Apply the selection method
                predictions, deltas, topTwo = self.apply_selectionMethod(seqAnnScores,selectionMethod, threshold_type)
                
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

            def apply_selectionMethod(self, seqAnnScores, selectionMethod,threshold_type):

                if selectionMethod == 'top_hit' and threshold_type is None:
                    def get_topTwo(row):
                        top2 = row.nlargest(2)
                        return pd.Series({
                            'keyValueOne': top2.iloc[0] if len(top2) > 0 else np.nan,
                            'keyValueOneHeader': top2.index[0] if len(top2) > 0 else np.nan,
                            'keyValueTwo': top2.iloc[1] if len(top2) > 1 else np.nan,
                            'keyValueTwoHeader': top2.index[1] if len(top2) > 1 else np.nan
                        })

                    keyValsDf = seqAnnScores.apply(get_topTwo, axis=1)
                    predictions = keyValsDf['keyValueOneHeader'].tolist()
                    deltas = (keyValsDf['keyValueOne'] - keyValsDf['keyValueTwo']).tolist()
                    topTwo = keyValsDf[['keyValueOne', 'keyValueTwo']]

                    return predictions, deltas, topTwo

                elif selectionMethod == 'top_hit' and threshold_type is not None:
                    def get_topTwo(row):
                        # Filter values based on threshold, keeping only values above the threshold
                        thresholds = row.index.map(self.threshold_dict).to_numpy()  # Map thresholds to each family (column) for this row
                        filtered_row = row.where(row >= thresholds, np.nan)

                        # Get the top two values and their column headers, handling cases with fewer than two valid entries
                        top2 = filtered_row.nlargest(2)
                        return pd.Series({
                            'keyValueOne': top2.iloc[0] if len(top2) > 0 else np.nan,
                            'keyValueOneHeader': top2.index[0] if len(top2) > 0 else np.nan,
                            'keyValueTwo': top2.iloc[1] if len(top2) > 1 else np.nan,
                            'keyValueTwoHeader': top2.index[1] if len(top2) > 1 else np.nan
                        })

                    # Apply the function to each row of seqAnnScores
                    keyValsDf = seqAnnScores.apply(get_topTwo, axis=1)

                    # Extracting the predictions and deltas for output
                    predictions = keyValsDf['keyValueOneHeader'].tolist()
                    deltas = (keyValsDf['keyValueOne'] - keyValsDf['keyValueTwo']).tolist()
                    topTwo = keyValsDf[['keyValueOne', 'keyValueTwo']]
                    return predictions, deltas, topTwo

                # The top two values compared here are: Greatest dist from threshold and Second greatest dist from Threshold
                elif selectionMethod == 'greatest_distance':
                    thresholds = seqAnnScores.columns.to_series().map(self.threshold_dict)
                    distances = seqAnnScores.subtract(thresholds, axis=1)
                    filtered_distances = distances.where(distances >= 0, np.nan)
                    partitioned_distances = np.partition(
                        np.nan_to_num(filtered_distances.values, nan=-np.inf), -2, axis=1
                    )
                    
                    topTwoIdx = np.argpartition(
                        np.nan_to_num(filtered_distances.values, nan=-np.inf), -2, axis=1
                    )[:, -2:]
                    
                    sortedTopTwoIdx = np.argsort(
                        filtered_distances.values[
                            np.arange(filtered_distances.shape[0])[:, None], topTwoIdx
                        ],
                        axis=1,
                    )[:, ::-1]
                    
                    topTwoDistances = np.take_along_axis(
                        filtered_distances.values,
                        topTwoIdx[np.arange(filtered_distances.shape[0])[:, None], sortedTopTwoIdx],
                        axis=1,
                    )
                    
                    topTwoScores = np.take_along_axis(
                        seqAnnScores.values,
                        topTwoIdx[np.arange(filtered_distances.shape[0])[:, None], sortedTopTwoIdx],
                        axis=1,
                    )
                    
                    topTwoHeaders = np.array(filtered_distances.columns)[
                        topTwoIdx[np.arange(filtered_distances.shape[0])[:, None], sortedTopTwoIdx]
                    ]
                    
                    flattened_headers = topTwoHeaders.flatten()
                    flattened_thresholds = thresholds[flattened_headers].values
                    topTwoThresholds = flattened_thresholds.reshape(topTwoHeaders.shape)
                    
                    keyValsDf = pd.DataFrame({
                        'keyValueOne': topTwoScores[:, 0],
                        'keyValueOneHeader': topTwoHeaders[:, 0],
                        'keyValueOneDistance': topTwoDistances[:, 0],
                        'keyValueOneThreshold': topTwoThresholds[:, 0],
                        'keyValueTwo': topTwoScores[:, 1],
                        'keyValueTwoHeader': topTwoHeaders[:, 1],
                        'keyValueTwoDistance': topTwoDistances[:, 1],
                        'keyValueTwoThreshold': topTwoThresholds[:, 1],
                    })
                    
                    deltas = (keyValsDf['keyValueOneDistance'] - keyValsDf['keyValueTwoDistance']).tolist()
                    predictions = keyValsDf['keyValueOneHeader'].tolist()
                    topTwo = keyValsDf[['keyValueOne', 'keyValueTwo']]
                    
                    return predictions, deltas, topTwo

        
                elif selectionMethod == 'combined_distance':
                    # Method 4: Combined method with dynamic delta calculation
                    # Set weights
                    weight_top = config["learnapp"].get("weight_top", 0.5)
                    weight_distance = config["learnapp"].get("weight_distance", 0.5)
                    thresholds = seqAnnScores.columns.to_series().map(self.threshold_dict)
                    distances = seqAnnScores - thresholds
                    positiveDistances = distances.where(distances >= 0, np.nan)
                    combinedScores = (seqAnnScores * weight_top) + (distances * weight_distance)
                    combinedScores = combinedScores.where(positiveDistances.notna(), np.nan)
                    topCombinedScores = combinedScores.max(axis=1)
                    predictions = combinedScores.idxmax(axis=1)
                    predictions = predictions.where(~topCombinedScores.isna(), None)

                    deltas = []
                    topTwoList = []

                    ## Method 1/2: Top hit with threshold
                    filteredScores = seqAnnScores.where(seqAnnScores >= thresholds, np.nan)
                    top_scores_method1 = filteredScores.max(axis=1)
                    top_families_method1 = filteredScores.idxmax(axis=1)

                    # Get second highest scores
                    temp_scores_method1 = filteredScores.apply(lambda row: row[row != row.max()], axis=1)
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
                            delta = top_scores_method1.iloc[idx] - secondScoresMethod1.iloc[idx]
                            # For topTwo, use top two scores from method 1
                            topTwo = [top_scores_method1.iloc[idx], secondScoresMethod1.iloc[idx]]
                        elif pred_family == topFamilyMethod3:
                            # Use delta from method 3: difference between top score and threshold
                            original_score = seqAnnScores.loc[seqAnnScores.index[idx], pred_family]
                            threshold = thresholds[pred_family]
                            delta = original_score - threshold
                            # For topTwo, use original score and threshold
                            topTwo = [original_score, threshold]
                        else:
                            # Use delta from combined scores: difference between top two combined scores
                            row_combinedScores = combinedScores.iloc[idx]
                            # Exclude the top prediction to find the second highest combined score
                            secondCombinedScore = row_combinedScores.drop(pred_family).max()
                            delta = topCombinedScores.iloc[idx] - secondCombinedScore
                            # For topTwo, use top two combined scores
                            topTwo = [topCombinedScores.iloc[idx], secondCombinedScore]

                        deltas.append(delta)
                        topTwoList.append(topTwo)

                    # Create a DataFrame for topTwo
                    topTwo = pd.DataFrame(topTwoList, columns=['keyValueOne', 'keyValueTwo'])

                    # Convert predictions to list
                    predictions = predictions.tolist()

                    return predictions, deltas, topTwo

                else:
                    raise ValueError(f"Invalid selection method: {selectionMethod}")


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
                rounded_deltas = [round(num, 2) if num is not None else 0 for num in deltas]

                # print(f"two_key_vals in generateDiffDataframe:\n {two_key_vals}")
                diffDataframe = pd.DataFrame({
                    "Top": two_key_vals['keyValueOne'],
                    "Second": two_key_vals['keyValueTwo'],
                    "Difference": rounded_deltas,
                    "Prediction": predictions,
                    "Actual": result,
                    "T/F": tf,
                    "Known/Unknown": known
                })

                diffDataframe['Difference'] = pd.to_numeric(diffDataframe['Difference'], errors='coerce')

                return diffDataframe

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
                    (valid_predictions["Known/Unknown"] == "Known") & (valid_predictions["T/F"] == "T")
                ]
                known_false_diffDataframe = valid_predictions[
                    (valid_predictions["Known/Unknown"] == "Known") & (valid_predictions["T/F"] == "F")
                ]

                trueCrosstab = pd.crosstab(
                    known_true_diffDataframe.Prediction, known_true_diffDataframe.Difference
                )
                falseCrosstab = pd.crosstab(
                    known_false_diffDataframe.Prediction, known_false_diffDataframe.Difference
                )
                return trueCrosstab, falseCrosstab

            def handleRunningCrosstabs(
                self, trueCrosstab, falseCrosstab, iteration
            ):
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
                self.trueRunningCrosstab, self.falseRunningCrosstab = self.trueRunningCrosstab.align(
                    self.falseRunningCrosstab, join='outer', axis=1, fill_value=0
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
                    trueCrosstab, falseCrosstab = self.createTrueFalseCrosstabs(diffDataframe)
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
                trueTotalDist = self.trueRunningCrosstab.sum(
                    numeric_only=True, axis=0
                )
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
                newIndex = pd.Index([round(i, 2) for i in pd.np.arange(0, 1.01, 0.01)], name="Difference")

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
                currentWeight = self.trueRunningCrosstab.values.sum() + self.falseRunningCrosstab.values.sum()

                # Prepare the updatedData DataFrame with the new columns
                updatedData = pd.DataFrame(newRatioDist, columns=["confidence"])
                updatedData["weight"] = currentWeight
                updatedData["total_sum"] = total_sum  # no interpolation, raw cumulative counts
                updatedData["inter_sum"] = inter_sum  # interpolated cumulative counts

                if self.confidenceData and len(self.confidenceData) == 1:
                    priorConf = pd.read_csv(self.confidenceData[0], index_col="Difference")
                    print(f"Prior Confidence Data:\n{priorConf}")

                    totalWeightPrior = priorConf["weight"]
                    kFactor = 1 + self.modifier * (currentWeight / (currentWeight + totalWeightPrior))
                    outWeight = totalWeightPrior + currentWeight
                    weightedCurrent = kFactor * currentWeight
                    totalWeight = totalWeightPrior + weightedCurrent
                    priorWeightedScore = priorConf["confidence"] * priorConf["weight"]
                    currentWeighted_score = updatedData["confidence"] * weightedCurrent

                    updatedConfidence = (priorWeightedScore + currentWeighted_score) / totalWeight
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
                        updatedData["cur_sum"] = priorConf["cur_sum"] + (self.trueRunningCrosstab.sum() + self.falseRunningCrosstab.sum())
                    else:
                        # If no prior sum, just create it now
                        updatedData["cur_sum"] = self.trueRunningCrosstab.sum() + self.falseRunningCrosstab.sum()

                    updatedData["weight"] = outWeight
                    print(f"Final Confidence Data\n{updatedData}")
                else:
                    print("Base confidence file not found or multiple files present. Only one file is allowed in baseConfidence.")
                    # If no merging, assign sum as needed
                    updatedData["cur_sum"] = self.trueRunningCrosstab.sum() + self.falseRunningCrosstab.sum()
                
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
                ratioDist, total_sum, inter_sum = self.computeRatioDistribution(trueDist, falseDist)
                globalConfidence = self.checkConfidenceMerge(ratioDist, total_sum, inter_sum)
                self.saveResults(
                    globalConfidence,
                    self.outputGlob
                )

        evaluator = Evaluator(
            input.evalApplyData,
            output.evalGlob,
            input.reverseDecoyStats,
            params.modifier,
            input.baseConfidence,
        )
        evaluator.executeAll()

        skm.utils.log_runtime(log[0], start_time)

rule copy_results_for_apply:
    input:
        kmer_counts_total = join(outDir, "learn", "kmer-counts-total.csv"),
        family_stats = join(outDir, "eval_conf", "family_summary_stats.csv"),
        global_conf_scores = join(outDir, "eval_conf", "global-confidence-scores.csv")
    output:
        kmer_counts_total=join(outDir, "apply_inputs","counts","kmer-counts-total.csv"),
        family_stats=join(outDir, "apply_inputs", "stats","family_summary_stats.csv"),
        global_conf_scores=join(outDir, "apply_inputs","confidence", "global-confidence-scores.csv")
    run:
        target_dir = os.path.join(outDir, "apply_inputs")
        os.makedirs(target_dir, exist_ok=True)
        shutil.copy(input.kmer_counts_total, output.kmer_counts_total)
        shutil.copy(input.family_stats, output.family_stats)
        shutil.copy(input.global_conf_scores, output.global_conf_scores)


        