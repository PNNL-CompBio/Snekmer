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
# reverseDecoy_files = glob(join("decoys", "*"))
# reverseDecoy_basenames = [skm.utils.split_file_ext(f)[0] for f in reverseDecoy_files]
# reverseDecoy_map = {
#     skm.utils.split_file_ext(f)[0]: skm.utils.split_file_ext(f)[1] for f in reverseDecoy_files
# }
# map extensions to basename (basename.ext.gz -> {basename: ext})
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

# # Add a new rule to reverse sequences
# rule reverse_sequences:
#     input:
#         fasta=lambda wildcards: join(
#             inputDir,
#             f"{wildcards.nb}.{faMap[wildcards.nb]}",
#         ),
#     output:
#         fasta_out=lambda wildcards: join("output","reversed", f"{wildcards.nb}_reversed.{faMap[wildcards.nb]}"),
#     run:
#         with open(input.fasta, 'r') as f_in, open(output.fasta_out, 'w') as f_out:
#             for record in SeqIO.parse(f_in, 'fasta'):
#                 reversed_seq = record.seq[::-1]  # Reverse the sequence
                
#                 # Create a new record
#                 reversed_record = record
#                 reversed_record.seq = reversed_seq
#                 # Modify the record id or description if needed
#                 reversed_record.id = record.id + '_reversed'
#                 reversed_record.description = record.description + ' reversed'
#                 SeqIO.write(reversed_record, f_out, 'fasta')

# define output files to be created by snekmer
rule all:
    input:
        expand(join("output", "vector", "{nb}.npz"), nb=FAS),
        # expand(join("output", "vector_decoy", "{dc}.npz"), dc=reverseDecoy_basenames),
        expand(join("output", "learn", "kmer-counts-{nb}.csv"), nb=FAS),
        join("output", "learn", "kmer-counts-total.csv"),
        expand(join("output", "fragmented", "{nb}.fasta"), nb=FAS)
        if config["learnapp"]["fragmentation"]
        else [],
        expand(join("output", "vector_frag", "nb}.npz"), nb=FAS)
        if config["learnapp"]["fragmentation"]
        else [],
        expand(
            join(
                "output",
                "eval_apply_sequences"
                if not config["learnapp"]["fragmentation"]
                else "eval_apply_frag",
                "seq-annotation-scores-{nb}.csv.gz",
            ),
            nb=FAS,
        ),
        # expand(
        #     join(
        #         "output",
        #         "eval_apply_reversed"
        #         if not config["learnapp"]["fragmentation"]
        #         else "eval_apply_frag",
        #         "seq-annotation-scores-{nb}.csv",
        #     ),
        #     nb=FAS,
        # ),
        # "output/evalConf/confidence-matrix.csv",
        # "output/evalConf/global-confidence-scores.csv",
        # "output/evalConf/max_scores_distribution.png",
        # "output/evalConf/max_scores_summary_stats.csv",
        # "output/evalConf/median_values_distribution.png",
        "output/eval_conf/family_summary_stats.csv",
        "output/eval_conf/global-confidence-scores.csv",


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
            fasta_out=join("output", "fragmented", "{nb}.fasta"),
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


# ruleorder:
#     vectorize > vectorize_decoy

use rule vectorize from kmerize with:
    input:
        fasta=lambda wildcards: join(
            "output" if wildcards.prefix == "vector_frag" else inputDir,
            "fragmented" if wildcards.prefix == "vector_frag" else "",
            f"{wildcards.nb}.{faMap[wildcards.nb]}",
        ),
    output:
        data=join("output", "{prefix}", "{nb}.npz"),
        kmerobj=join("output", "kmerize_{prefix}", "{nb}.kmers"),
    log:
        join("output", "{prefix}_kmerize", "log", "{nb}.log"),

# use rule vectorize from kmerize as vectorize_decoy with:
#     input:
#         # fasta=lambda wildcards: expand(join("decoys", "{dc}.fa"), dc=reverseDecoy_basenames)[0],
#         fasta=lambda wildcards: join(
#             "decoys",
#             f"{wildcards.dc}.{reverseDecoy_map[wildcards.dc]}",
#         ),
#     output:
#         data=join("output", "vector_decoy", "{dc}.npz"),  # Output path updated
#         kmerobj=join("output", "kmerize_decoy", "{dc}.kmers"),  # Output path updated
#     log:
#         join("output", "vector_decoy", "log", "{dc}.log"),  # Log path updated


# WORKFLOW to learn kmer associations
# collect all seq files and generate mega-cluster
rule learn:
    input:
        data="output/vector/{nb}.npz",
        annotation=expand("{an}", an=annotFiles),
    output:
        counts="output/learn/kmer-counts-{nb}.csv",
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
                anns = annotations["Family"].tolist()
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
                    x = re.findall(r"\|(.*?)\|", seqid)[0]
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
                out_name = "output/learn/kmer-counts-" + str(inputData)[14:-4] + ".csv"
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
        counts=expand(join("output", "learn", "kmer-counts-{nb}.csv"), nb=FAS),
        baseCounts=expand("{bf}", bf=baseCounts),
    output:
        totals=join("output", "learn", "kmer-counts-total.csv"),
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
            "output",
            "vector"
            if config["learnapp"]["fragmentation"] == False
            else "vector_frag",
            "{nb}.npz",
        ),
        annotation=expand("{an}", an=annotFiles),
        compareAssociations=join("output", "learn", "kmer-counts-total.csv"),

    output:
        apply=join(
            "output",
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
                finalMatrixWithScores = finalMatrixWithScores.round(3)  #Adding rounding for storage saving.
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
                # finalMatrixWithScores = self.filterTopTwoValues(finalMatrixWithScores)  # We can't do this because we don't know the thresholds at this point yet.
                self.writeOutput(finalMatrixWithScores)
                # Add next things here
                

        analysis = KmerCompare(
            input.compareAssociations, input.annotation, input.data, output.apply
        )
        analysis.executeAll(config)
        skm.utils.log_runtime(log[0], start_time)


rule reverseDecoy_evaluations:
    input:
        evalApplyData=expand(
            join(
                "output",
                "eval_apply_reversed"
                if config["learnapp"]["fragmentation"] == False
                else "eval_apply_frag",
                "seq-annotation-scores-{nb}.csv.gz",
            ),
            nb=FAS,
        )
    output:
        familyStats="output/eval_conf/family_summary_stats.csv",
    run:

        # def collectFamilyStatistics(filename, existingStats=None):
        #     """
        #     Reads a CSV file in chunks and updates statistics for each family (column).

        #     Args:
        #         filename (str): Path to the CSV file.
        #         existingStats (dict): Existing statistics to update.

        #     Returns:
        #         dict: Updated statistics for each family.
        #     """
        #     chunk_size = 10000  # Adjust based on available memory
        #     if existingStats is None:
        #         existingStats = {}

        #     # Initialize variables to compute running statistics
        #     for chunk in pd.read_csv(filename, chunksize=chunk_size, engine="c"):
        #         # Exclude the last column if it's the sequence name
        #         families = chunk.columns[:-1]

        #         # Process each family (column)
        #         for family in families:
        #             values = chunk[family].dropna().astype(float)

        #             if family not in existingStats:
        #                 existingStats[family] = {
        #                     'count': 0,
        #                     'mean': 0.0,
        #                     'M2': 0.0,  # Sum of squares of differences from the current mean
        #                     'min': np.inf,
        #                     'max': -np.inf,
        #                     'values_for_percentiles': []
        #                 }

        #             stats = existingStats[family]

        #             # Update count
        #             n = len(values)
        #             if n == 0:
        #                 continue
        #             stats['count'] += n

        #             # Update mean and M2 for standard deviation (Welford's algorithm)
        #             delta = values - stats['mean']
        #             stats['mean'] += delta.sum() / stats['count']
        #             delta2 = values - stats['mean']
        #             stats['M2'] += (delta * delta2).sum()

        #             # Update min and max
        #             stats['min'] = min(stats['min'], values.min())
        #             stats['max'] = max(stats['max'], values.max())

        #             # For percentiles, we can store a sample if data is too large
        #             stats['values_for_percentiles'].extend(values.tolist())
        #             # Optionally limit the size of the list to save memory
        #             if len(stats['values_for_percentiles']) > 100000:
        #                 stats['values_for_percentiles'] = np.random.choice(
        #                     stats['values_for_percentiles'], 100000, replace=False
        #                 ).tolist()

        #         # Clean up to free memory
        #         del chunk
                
        #     return existingStats

        # def generateFamilyStatistics(combinedStats):
        #     """
        #     Generates statistics for each family using the combined statistics.

        #     Args:
        #         combinedStats (dict): Combined statistics for each family.

        #     Returns:
        #         pd.DataFrame: DataFrame containing the statistics for each family.
        #     """
        #     statsData = {
        #         'Family': [],
        #         'Mean': [],
        #         'Std Dev': [],
        #         'Min': [],
        #         '10th Percentile': [],
        #         '20th Percentile': [],
        #         '25th Percentile': [],
        #         '30th Percentile': [],
        #         '40th Percentile': [],
        #         'Median': [],
        #         '60th Percentile': [],
        #         '70th Percentile': [],
        #         '75th Percentile': [],
        #         '80th Percentile': [],
        #         '90th Percentile': [],
        #         'Max': [],
        #         '1 Std Dev Above': [],
        #         '1 Std Dev Below': [],
        #         '2 Std Dev Above': [],
        #         '2 Std Dev Below': [],
        #     }

        #     for family, stats in combinedStats.items():
        #         count = stats['count']
        #         mean = stats['mean']
        #         std_dev = np.sqrt(stats['M2'] / (count - 1)) if count > 1 else 0.0

        #         # Use the stored values to compute percentiles
        #         values = np.array(stats['values_for_percentiles'])
        #         percentiles = np.percentile(values, [10, 20, 25, 30, 40, 50, 60, 70, 75, 80, 90])

        #         statsData['Family'].append(family)
        #         statsData['Mean'].append(round(mean, 3))
        #         statsData['Std Dev'].append(round(std_dev, 3))
        #         statsData['Min'].append(round(stats['min'], 3))
        #         statsData['10th Percentile'].append(round(percentiles[0], 3))
        #         statsData['20th Percentile'].append(round(percentiles[1], 3))
        #         statsData['25th Percentile'].append(round(percentiles[2], 3))
        #         statsData['30th Percentile'].append(round(percentiles[3], 3))
        #         statsData['40th Percentile'].append(round(percentiles[4], 3))
        #         statsData['Median'].append(round(percentiles[5], 3))
        #         statsData['60th Percentile'].append(round(percentiles[6], 3))
        #         statsData['70th Percentile'].append(round(percentiles[7], 3))
        #         statsData['75th Percentile'].append(round(percentiles[8], 3))
        #         statsData['80th Percentile'].append(round(percentiles[9], 3))
        #         statsData['90th Percentile'].append(round(percentiles[10], 3))
        #         statsData['Max'].append(round(stats['max'], 3))
        #         statsData['1 Std Dev Above'].append(round(mean + std_dev, 3))
        #         statsData['1 Std Dev Below'].append(round(mean - std_dev, 3))
        #         statsData['2 Std Dev Above'].append(round(mean + 2 * std_dev, 3))
        #         statsData['2 Std Dev Below'].append(round(mean - 2 * std_dev, 3))

        #     return pd.DataFrame(statsData)


        # def collectFamilyStatistics(filename, existingStats=None):
        #     """
        #     Reads a CSV file in chunks and updates statistics for each family (column).

        #     Args:
        #         filename (str): Path to the CSV file.
        #         existingStats (dict): Existing statistics to update.

        #     Returns:
        #         dict: Updated statistics for each family.
        #     """
        #     chunk_size = 10000  # Adjust based on available memory
        #     if existingStats is None:
        #         existingStats = {}

        #     for chunk in pd.read_csv(filename, chunksize=chunk_size, engine="c"):
        #         families = chunk.columns[:-1]  # Exclude the last column if it's the sequence name

        #         # Process each family (column)
        #         for family in families:
        #             values = chunk[family].dropna().astype(float).values
        #             if family not in existingStats:
        #                 existingStats[family] = {
        #                     'count': 0,
        #                     'sum': 0.0,
        #                     'sumSqr': 0.0,
        #                     'min': np.inf,
        #                     'max': -np.inf,
        #                     'values_for_percentiles': []
        #                 }

        #             stats = existingStats[family]

        #             n = len(values)
        #             if n == 0:
        #                 continue
        #             stats['count'] += n
        #             stats['sum'] += values.sum()
        #             stats['sumSqr'] += (values ** 2).sum()

        #             stats['min'] = min(stats['min'], values.min())
        #             stats['max'] = max(stats['max'], values.max())  

        #             # For percentiles, we need to use sampling to reduce memory.
        #             stats['values_for_percentiles'].extend(values.tolist())
        #             if len(stats['values_for_percentiles']) > 100000:
        #                 stats['values_for_percentiles'] = np.random.choice(
        #                     stats['values_for_percentiles'], 100000, replace=False
        #                 ).tolist()

        #         del chunk  # Clean up to free memory

        #     return existingStats



        def collectFamilyStatistics(filename, existingStats=None):
            """
            Reads a CSV file in chunks and updates statistics for each family (column).

            Args:
                filename (str): Path to the CSV file.
                existingStats (dict): Existing statistics to update.

            Returns:
                dict: Updated statistics for each family.
            """
            chunk_size = 10000  # Adjust based on available memory
            reservoir_size = 100000  # Size of the reservoir for percentiles

            if existingStats is None:
                existingStats = {}

            for chunk in pd.read_csv(filename, chunksize=chunk_size, engine="c"):
                families = chunk.columns[:-1]  # Exclude the last column if it's the sequence name

                # Process each family (column)
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

                del chunk  # Clean up to free memory

            return existingStats

        def generateFamilyStatistics(combinedStats):
            """
            Generates statistics for each family using the combined statistics.

            Args:
                combinedStats (dict): Combined statistics for each family.

            Returns:
                pd.DataFrame: DataFrame containing the statistics for each family.
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
                mean = sum_ / n
                variance = (sumSqr - (sum_ ** 2) / n) / (n - 1) if n > 1 else 0.0
                std_dev = np.sqrt(variance)

                # Use all stored values to compute percentiles
                values = np.array(stats['values_for_percentiles'])
                percentiles = np.percentile(values, [10, 20, 25, 30, 40, 50, 60, 70, 75, 80, 90])

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


        # Main execution
        # Step 1: Collect statistics for each family from all input files
        
        combinedStats = {}

        for filename in input.evalApplyData:
            combinedStats = collectFamilyStatistics(filename, existingStats=combinedStats)

        # Step 2: Generate summary statistics for each family
        familyStatisticsDf = generateFamilyStatistics(combinedStats)

        # Save family statistics to CSV
        familyStatisticsDf.to_csv(output.familyStats, index=False)

        #Needs lots of memory
        # def collect_family_values(filename):
        #     """
        #     Reads a CSV file and returns a dictionary with all values for each family (column).
        #     """
        #     # Read the CSV file using pandas
        #     seqAnnScores = pd.read_csv(
        #         filename,
        #         index_col=None,
        #         header=0,
        #         engine="c"
        #     )

        #     # Dictionary to store values for each family (column)
        #     family_values = {family: [] for family in seqAnnScores.columns[:-1]}  # Exclude the last column (sequence name)

        #     # Iterate over the rows in the DataFrame
        #     for index, row in seqAnnScores.iterrows():
        #         # Iterate over all families (columns except the last one)
        #         for family in family_values.keys():
        #             value = row[family]
        #             if not pd.isna(value):  # Check if the value is not NaN
        #                 family_values[family].append(float(value))

        #     return family_values

        # def generateFamilyStatistics(family_values):
        #     """
        #     Generates statistics (max, min, mean, median, percentiles, and standard deviations) for each family, rounded to 3 decimal places.
        #     """
        #     statsData = {
        #         'Family': [],
        #         'Mean': [],
        #         '1 Std Dev Above': [],
        #         '1 Std Dev Below': [],
        #         '2 Std Dev Above': [],
        #         '2 Std Dev Below': [],
        #         'Min': [],
        #         '10th Percentile': [],
        #         '20th Percentile': [],
        #         '25th Percentile': [],
        #         '30th Percentile': [],
        #         '40th Percentile': [],
        #         'Median': [],
        #         '60th Percentile': [],
        #         '70th Percentile': [],
        #         '75th Percentile': [],
        #         '80th Percentile': [],
        #         '90th Percentile': [],
        #         'Max': [],
        #     }

        #     for family, values in family_values.items():
        #         mean_val = np.mean(values)
        #         std_dev = np.std(values)

        #         statsData['Family'].append(family)
        #         statsData['Mean'].append(round(mean_val, 3))
        #         statsData['Min'].append(round(np.min(values), 3))
        #         statsData['10th Percentile'].append(round(np.percentile(values, 10), 3))
        #         statsData['20th Percentile'].append(round(np.percentile(values, 20), 3))
        #         statsData['25th Percentile'].append(round(np.percentile(values, 25), 3))
        #         statsData['30th Percentile'].append(round(np.percentile(values, 30), 3))
        #         statsData['40th Percentile'].append(round(np.percentile(values, 40), 3))
        #         statsData['Median'].append(round(np.median(values), 3))
        #         statsData['60th Percentile'].append(round(np.percentile(values, 60), 3))
        #         statsData['70th Percentile'].append(round(np.percentile(values, 70), 3))
        #         statsData['75th Percentile'].append(round(np.percentile(values, 75), 3))
        #         statsData['80th Percentile'].append(round(np.percentile(values, 80), 3))
        #         statsData['90th Percentile'].append(round(np.percentile(values, 90), 3))
        #         statsData['Max'].append(round(np.max(values), 3))
        #         statsData['1 Std Dev Above'].append(round(mean_val + std_dev, 3))
        #         statsData['1 Std Dev Below'].append(round(mean_val - std_dev, 3))
        #         statsData['2 Std Dev Above'].append(round(mean_val + 2 * std_dev, 3))
        #         statsData['2 Std Dev Below'].append(round(mean_val - 2 * std_dev, 3))

        #     return pd.DataFrame(statsData)

        # def plot_distribution(data, column_name, output_file):
        #     """
        #     Plots the distribution of a specified column in the given data.
        #     """
        #     plt.figure(figsize=(10, 6))
        #     plt.hist(data[column_name], bins=30, alpha=0.75, edgecolor='black')
        #     plt.title(f'Distribution of {column_name} Values for Each Family')
        #     plt.xlabel(f'{column_name} Value')
        #     plt.ylabel('Frequency')
        #     plt.grid(axis='y', alpha=0.75)
        #     plt.savefig(output_file)
        #     plt.close()

        # # Step 1: Collect values for each family from all input files
        # combined_family_values = {}

        # for filename in input.evalApplyData:
        #     family_values = collect_family_values(filename)

        #     # Combine values from multiple files if needed
        #     for family, values in family_values.items():
        #         if family not in combined_family_values:
        #             combined_family_values[family] = []
        #         combined_family_values[family].extend(values)

        # # Step 2: Generate summary statistics for each family
        # familyStatisticsDf = generateFamilyStatistics(combined_family_values)

        # # Save family statistics to CSV
        # familyStatisticsDf.to_csv(output.familyStats, index=False)




###
###
###

# This is working - but way too slow.


# rule combined_eval_apply_reverse_seqs_and_evaluations:
#     input:
#         data=expand(
#             join(
#                 "output",
#                 "vector"
#                 if config["learnapp"]["fragmentation"] == False
#                 else "vector_frag",
#                 "{nb}.npz",
#             ),
#             nb=FAS,
#         ),
#         annotation=annotFiles,
#         compareAssociations=join("output", "learn", "kmer-counts-total.csv"),
#     output:
#         familyStats="output/eval_conf/family_summary_stats.csv",
#     log:
#         "output/eval_apply_reversed/log/combined.log",
#     run:
#         from datetime import datetime
#         import pandas as pd
#         import numpy as np
#         import sklearn.metrics.pairwise
#         import sys
#         import itertools
#         import os

#         import snekmer as skm

#         start_time = datetime.now()
#         with open(log[0], "a") as f:
#             f.write(f"start time:\t{start_time}\n")

#         class KmerCompare:
#             """
#             Initializes the KmerCompare object.

#             This object is designed to compare kmer counts with provided annotations.

#             Attributes:
#                 compareAssociations (str): Path to a CSV file containing kmer counts totals matrix.
#                 annotationFiles (list): List of paths to files containing sequence annotations.
#                 inputData (list): List of paths to input data for kmer analysis.
#                 seqAnnot (dict): Dictionary mapping sequence IDs to annotations.
#                 kmerCountTotals (DataFrame): DataFrame of kmer count totals from compareAssociations.
#                 kmerList (list): List of unique kmers found in input data.
#                 combinedStats (dict): Dictionary to store family statistics.
#             """

#             def __init__(self, compareAssociations, annotationFiles, inputData):
#                 self.compareAssociations = compareAssociations
#                 self.annotationFiles = annotationFiles
#                 self.inputData = inputData
#                 self.seqAnnot = {}
#                 self.kmerCountTotals = None
#                 self.kmerList = []
#                 self.combinedStats = {}

#             def generateInputs(self):
#                 """
#                 Generates the necessary inputs for comparison.

#                 Loads kmer counts and annotations into appropriate data structures.
#                 """
#                 self.kmerCountTotals = pd.read_csv(
#                     str(self.compareAssociations),
#                     index_col="__index_level_0__",
#                     header=0,
#                     engine="c",
#                 )
#                 for f in self.annotationFiles:
#                     annotation_df = pd.read_table(f)
#                     seqs = annotation_df["id"].tolist()
#                     anns = annotation_df["Family"].tolist()
#                     for seqid, ann in zip(seqs, anns):
#                         self.seqAnnot[seqid] = ann

#             def executeAll(self, config):
#                 """
#                 Executes all the comparison steps in sequence.

#                 This includes:
#                     1. Loading inputs.
#                     2. Processing each .npz file.
#                     3. Generating kmer counts for reversed sequences.
#                     4. Calculating cosine similarity and updating family statistics.
#                     5. Generating and writing family statistics.
#                 """
#                 self.generateInputs()

#                 for npz_file in self.inputData:
#                     kmerList, df = skm.io.load_npz(npz_file)
#                     self.kmerList = kmerList[0]
#                     seqids = df["sequence_id"]
#                     sequences = df["sequence"]
#                     k_len = len(self.kmerList[0])

#                     for seq_id, seq in zip(seqids, sequences):
#                         v = seq[::-1]

#                         kCounts = {}
#                         items = [
#                             v[i : i + k_len]
#                             for i in range(len(v) - k_len + 1)
#                         ]
#                         for kmer in items:
#                             kCounts[kmer] = kCounts.get(kmer, 0) + 1
#                         store = [kCounts.get(kmer, 0) for kmer in self.kmerList]

#                         # Compute cosine similarity for this sequence
#                         kmerCounts = pd.DataFrame(
#                             [store], index=[seq_id], columns=self.kmerList
#                         )
#                         kmerCounts = self.matchkmerCountsFormat(kmerCounts)
#                         cosine_scores = self.calculateCosineSimilarity(kmerCounts)
#                         self.updateFamilyStatistics(cosine_scores)

#                 familyStatisticsDf = self.generateFamilyStatistics()
#                 familyStatisticsDf.to_csv(output.familyStats, index=False)

#             def matchkmerCountsFormat(self, kmerCounts):
#                 """
#                 Matches the format of the provided kmer counts DataFrame to the format of the
#                 comparison data. Ensures columns align correctly.
#                 """
#                 kmerCounts = kmerCounts.reindex(
#                     columns=self.kmerCountTotals.columns, fill_value=0
#                 )
#                 return kmerCounts

#             def calculateCosineSimilarity(self, kmerCounts):
#                 """
#                 Calculates the cosine similarity between the input kmer counts and comparison data.

#                 Returns:
#                     Series: A Series of cosine similarity scores.
#                 """
#                 cosine_scores = sklearn.metrics.pairwise.cosine_similarity(
#                     self.kmerCountTotals, kmerCounts
#                 ).flatten()
#                 cosine_scores_series = pd.Series(
#                     cosine_scores, index=self.kmerCountTotals.index
#                 )
#                 return cosine_scores_series

#             def updateFamilyStatistics(self, cosine_scores):
#                 """
#                 Updates the family statistics with the provided cosine similarity scores.

#                 Args:
#                     cosine_scores (Series): Cosine similarity scores for one sequence.
#                 """
#                 for family, score in cosine_scores.items():
#                     if family not in self.combinedStats:
#                         self.combinedStats[family] = {
#                             "count": 0,
#                             "sum": 0.0,
#                             "sumSqr": 0.0,
#                             "min": np.inf,
#                             "max": -np.inf,
#                             "values_for_percentiles": [],
#                         }
#                     stats = self.combinedStats[family]
#                     stats["count"] += 1
#                     stats["sum"] += score
#                     stats["sumSqr"] += score ** 2
#                     stats["min"] = min(stats["min"], score)
#                     stats["max"] = max(stats["max"], score)
#                     stats["values_for_percentiles"].append(score)
#                     # Optionally limit the size of values_for_percentiles to save memory
#                     if len(stats["values_for_percentiles"]) > 100000:
#                         stats["values_for_percentiles"] = np.random.choice(
#                             stats["values_for_percentiles"], 100000, replace=False
#                         ).tolist()

#             def generateFamilyStatistics(self):
#                 """
#                 Generates statistics for each family using the combined statistics.

#                 Returns:
#                     DataFrame: DataFrame containing the statistics for each family.
#                 """
#                 statsData = {
#                     "Family": [],
#                     "Mean": [],
#                     "Std Dev": [],
#                     "Min": [],
#                     "10th Percentile": [],
#                     "20th Percentile": [],
#                     "25th Percentile": [],
#                     "30th Percentile": [],
#                     "40th Percentile": [],
#                     "Median": [],
#                     "60th Percentile": [],
#                     "70th Percentile": [],
#                     "75th Percentile": [],
#                     "80th Percentile": [],
#                     "90th Percentile": [],
#                     "Max": [],
#                     "1 Std Dev Above": [],
#                     "1 Std Dev Below": [],
#                     "2 Std Dev Above": [],
#                     "2 Std Dev Below": [],
#                 }
#                 for family, stats in self.combinedStats.items():
#                     n = stats["count"]
#                     sum_ = stats["sum"]
#                     sumSqr = stats["sumSqr"]
#                     mean = sum_ / n
#                     variance = (
#                         (sumSqr - (sum_ ** 2) / n) / (n - 1) if n > 1 else 0.0
#                     )
#                     std_dev = np.sqrt(variance)
#                     values = np.array(stats["values_for_percentiles"])
#                     percentiles = np.percentile(
#                         values, [10, 20, 25, 30, 40, 50, 60, 70, 75, 80, 90]
#                     )
#                     statsData["Family"].append(family)
#                     statsData["Mean"].append(round(mean, 3))
#                     statsData["Std Dev"].append(round(std_dev, 3))
#                     statsData["Min"].append(round(stats["min"], 3))
#                     statsData["10th Percentile"].append(round(percentiles[0], 3))
#                     statsData["20th Percentile"].append(round(percentiles[1], 3))
#                     statsData["25th Percentile"].append(round(percentiles[2], 3))
#                     statsData["30th Percentile"].append(round(percentiles[3], 3))
#                     statsData["40th Percentile"].append(round(percentiles[4], 3))
#                     statsData["Median"].append(round(percentiles[5], 3))
#                     statsData["60th Percentile"].append(round(percentiles[6], 3))
#                     statsData["70th Percentile"].append(round(percentiles[7], 3))
#                     statsData["75th Percentile"].append(round(percentiles[8], 3))
#                     statsData["80th Percentile"].append(round(percentiles[9], 3))
#                     statsData["90th Percentile"].append(round(percentiles[10], 3))
#                     statsData["Max"].append(round(stats["max"], 3))
#                     statsData["1 Std Dev Above"].append(round(mean + std_dev, 3))
#                     statsData["1 Std Dev Below"].append(round(mean - std_dev, 3))
#                     statsData["2 Std Dev Above"].append(round(mean + 2 * std_dev, 3))
#                     statsData["2 Std Dev Below"].append(round(mean - 2 * std_dev, 3))
#                 familyStatisticsDf = pd.DataFrame(statsData)
#                 return familyStatisticsDf

#         # Convert inputs to standard Python types
#         input_data = [str(f) for f in input.data]
#         annotation_files = [str(f) for f in input.annotation]
#         compare_associations = str(input.compareAssociations)

#         analysis = KmerCompare(
#             compare_associations, annotation_files, input_data
#         )
#         analysis.executeAll(config)
#         # Log runtime
#         end_time = datetime.now()
#         runtime = end_time - start_time
#         with open(log[0], "a") as f:
#             f.write(f"end time:\t{end_time}\n")
#             f.write(f"runtime:\t{runtime}\n")

###
###
###
# Version 3

# rule eval_apply_reverse_seqs:
#     input:
#         data=join(
#             "output",
#             "vector"
#             if config["learnapp"]["fragmentation"] == False
#             else "vector_frag",
#             "{nb}.npz",
#         ),
#         annotation=expand("{an}", an=annotFiles),
#         compareAssociations=join("output", "learn", "kmer-counts-total.csv"),
#     output:
#         chunkStats=join(
#             "output",
#             "eval_apply_reversed"
#             if config["learnapp"]["fragmentation"] == False
#             else "eval_apply_frag",
#             "chunk-stats-{nb}.pkl",
#         ),
#     log:
#         join(outDir, "eval_apply_reversed", "log", "{nb}.log"),
#     run:
#         import pickle
#         from datetime import datetime
#         import pandas as pd
#         import numpy as np
#         import sklearn.metrics.pairwise
#         import sys
#         import itertools
#         import os
#         import snekmer as skm

#         start_time = datetime.now()
#         with open(log[0], "a") as f:
#             f.write(f"start time:\t{start_time}\n")

#         class KmerCompare:
#             """
#             Initializes the KmerCompare object.

#             This object is designed to compare kmer counts with provided annotations.

#             Attributes:
#                 compareAssociations (str): Path to a CSV file containing kmer counts totals matrix.
#                 annotationFiles (list): List of paths to files containing sequence annotations.
#                 inputData (str): Path to input data for kmer analysis.
#                 outputPath (str): Path to save the result.
#                 seqAnnot (dict): Dictionary mapping sequence IDs to annotations.
#                 kmerCountTotals (DataFrame): DataFrame of kmer count totals from compareAssociations.
#                 kmerList (list): List of unique kmers found in input data.
#             """

#             def __init__(
#                 self, compareAssociations, annotationFiles, inputData, outputPath
#             ):
#                 self.compareAssociations = compareAssociations
#                 self.annotationFiles = annotationFiles
#                 self.inputData = inputData
#                 self.outputPath = outputPath
#                 self.seqAnnot = {}
#                 self.kmerCountTotals = None
#                 self.kmerList = []
#                 self.kmerTotals = []

#             def generateInputs(self):
#                 """
#                 Generates the necessary inputs for comparison.

#                 Loads kmer counts and annotations into appropriate data structures.
#                 """
#                 self.kmerCountTotals = pd.read_csv(
#                     str(self.compareAssociations),
#                     index_col="__index_level_0__",
#                     header=0,
#                     engine="c",
#                 )
#                 for f in self.annotationFiles:
#                     annotation_df = pd.read_table(f)
#                     seqs = annotation_df["id"].tolist()
#                     anns = annotation_df["Family"].tolist()
#                     for seqid, ann in zip(seqs, anns):
#                         self.seqAnnot[seqid] = ann

#             def executeAll(self, config):
#                 """
#                 Executes all the comparison steps in sequence.

#                 This includes:
#                     1. Loading inputs.
#                     2. Processing each sequence in the input data.
#                     3. Generating kmer counts for reversed sequences.
#                     4. Calculating cosine similarity and updating family statistics.
#                     5. Writing per-chunk family statistics.
#                 """
#                 self.generateInputs()
#                 # Initialize chunkStats
#                 chunkStats = {}
#                 # Load the input data
#                 kmerList, df = skm.io.load_npz(self.inputData)
#                 self.kmerList = kmerList[0]
#                 seqids = df["sequence_id"]
#                 sequences = df["sequence"]
#                 k_len = len(self.kmerList[0])

#                 for seq_id, seq in zip(seqids, sequences):
#                     v = seq[::-1]
#                     # Generate kmer counts
#                     kCounts = {}
#                     items = [v[i : i + k_len] for i in range(len(v) - k_len + 1)]
#                     for kmer in items:
#                         kCounts[kmer] = kCounts.get(kmer, 0) + 1
#                     # Create kmerCounts dataframe
#                     kmerCounts = pd.DataFrame([kCounts], index=[seq_id])
#                     kmerCounts = kmerCounts.reindex(columns=self.kmerList, fill_value=0)
#                     # Match format
#                     kmerCounts = self.matchkmerCountsFormat(kmerCounts)
#                     # Compute cosine similarities
#                     cosine_scores = self.calculateCosineSimilarity(kmerCounts)
#                     # Update chunkStats
#                     for family, score in cosine_scores.items():
#                         if family not in chunkStats:
#                             chunkStats[family] = {
#                                 'count': 0,
#                                 'sum': 0.0,
#                                 'sumSqr': 0.0,
#                                 'min': np.inf,
#                                 'max': -np.inf,
#                                 'values_for_percentiles': []
#                             }
#                         stats = chunkStats[family]
#                         stats['count'] += 1
#                         stats['sum'] += score
#                         stats['sumSqr'] += score ** 2
#                         stats['min'] = min(stats['min'], score)
#                         stats['max'] = max(stats['max'], score)
#                         stats['values_for_percentiles'].append(score)
#                         if len(stats['values_for_percentiles']) > 10000:
#                             stats['values_for_percentiles'] = np.random.choice(
#                                 stats['values_for_percentiles'], 10000, replace=False
#                             ).tolist()
#                 # Write chunkStats to output file
#                 self.writeChunkStats(chunkStats)

#             def matchkmerCountsFormat(self, kmerCounts):
#                 """
#                 Matches the format of the provided kmer counts DataFrame to the format of the
#                 comparison data. Ensures columns align correctly.
#                 """
#                 kmerCounts = kmerCounts.reindex(
#                     columns=self.kmerCountTotals.columns, fill_value=0
#                 )
#                 return kmerCounts

#             def calculateCosineSimilarity(self, kmerCounts):
#                 """
#                 Calculates the cosine similarity between the input kmer counts and comparison data.

#                 Returns:
#                     Series: A Series of cosine similarity scores.
#                 """
#                 cosine_scores = sklearn.metrics.pairwise.cosine_similarity(
#                     self.kmerCountTotals, kmerCounts
#                 ).flatten()
#                 cosine_scores_series = pd.Series(
#                     cosine_scores, index=self.kmerCountTotals.index
#                 )
#                 return cosine_scores_series

#             def writeChunkStats(self, chunkStats):
#                 """
#                 Writes the per-chunk family statistics to a pickle file.
#                 """
#                 with open(self.outputPath, 'wb') as f:
#                     pickle.dump(chunkStats, f)

#         analysis = KmerCompare(
#             input.compareAssociations, input.annotation, input.data, output.chunkStats
#         )
#         analysis.executeAll(config)
#         skm.utils.log_runtime(log[0], start_time)

# rule reverseDecoy_evaluations:
#     input:
#         chunkStatsData=expand(
#             join(
#                 "output",
#                 "eval_apply_reversed"
#                 if config["learnapp"]["fragmentation"] == False
#                 else "eval_apply_frag",
#                 "chunk-stats-{nb}.pkl",
#             ),
#             nb=FAS,
#         ),
#     output:
#         familyStats="output/eval_conf/family_summary_stats.csv",
#     run:
#         import pickle
#         import pandas as pd
#         import numpy as np

#         def combineChunkStats(chunkStatsFiles):
#             combinedStats = {}
#             for filename in chunkStatsFiles:
#                 with open(filename, 'rb') as f:
#                     chunkStats = pickle.load(f)
#                 for family, stats in chunkStats.items():
#                     if family not in combinedStats:
#                         combinedStats[family] = {
#                             'count': 0,
#                             'sum': 0.0,
#                             'sumSqr': 0.0,
#                             'min': np.inf,
#                             'max': -np.inf,
#                             'values_for_percentiles': []
#                         }
#                     combined = combinedStats[family]
#                     combined['count'] += stats['count']
#                     combined['sum'] += stats['sum']
#                     combined['sumSqr'] += stats['sumSqr']
#                     combined['min'] = min(combined['min'], stats['min'])
#                     combined['max'] = max(combined['max'], stats['max'])
#                     combined['values_for_percentiles'].extend(stats['values_for_percentiles'])
#                     if len(combined['values_for_percentiles']) > 100000:
#                         combined['values_for_percentiles'] = np.random.choice(
#                             combined['values_for_percentiles'], 100000, replace=False
#                         ).tolist()
#             return combinedStats

#         def generateFamilyStatistics(combinedStats):
#             """
#             Generates statistics for each family using the combined statistics.

#             Args:
#                 combinedStats (dict): Combined statistics for each family.

#             Returns:
#                 pd.DataFrame: DataFrame containing the statistics for each family.
#             """
#             statsData = {
#                 'Family': [],
#                 'Mean': [],
#                 'Std Dev': [],
#                 'Min': [],
#                 '10th Percentile': [],
#                 '20th Percentile': [],
#                 '25th Percentile': [],
#                 '30th Percentile': [],
#                 '40th Percentile': [],
#                 'Median': [],
#                 '60th Percentile': [],
#                 '70th Percentile': [],
#                 '75th Percentile': [],
#                 '80th Percentile': [],
#                 '90th Percentile': [],
#                 'Max': [],
#                 '1 Std Dev Above': [],
#                 '1 Std Dev Below': [],
#                 '2 Std Dev Above': [],
#                 '2 Std Dev Below': [],
#             }

#             for family, stats in combinedStats.items():
#                 n = stats['count']
#                 sum_ = stats['sum']
#                 sumSqr = stats['sumSqr']
#                 mean = sum_ / n
#                 variance = (sumSqr - (sum_ ** 2) / n) / (n - 1) if n > 1 else 0.0
#                 std_dev = np.sqrt(variance)

#                 # Use all stored values to compute percentiles
#                 values = np.array(stats['values_for_percentiles'])
#                 percentiles = np.percentile(values, [10, 20, 25, 30, 40, 50, 60, 70, 75, 80, 90])

#                 statsData['Family'].append(family)
#                 statsData['Mean'].append(round(mean, 3))
#                 statsData['Std Dev'].append(round(std_dev, 3))
#                 statsData['Min'].append(round(stats['min'], 3))
#                 statsData['10th Percentile'].append(round(percentiles[0], 3))
#                 statsData['20th Percentile'].append(round(percentiles[1], 3))
#                 statsData['25th Percentile'].append(round(percentiles[2], 3))
#                 statsData['30th Percentile'].append(round(percentiles[3], 3))
#                 statsData['40th Percentile'].append(round(percentiles[4], 3))
#                 statsData['Median'].append(round(percentiles[5], 3))
#                 statsData['60th Percentile'].append(round(percentiles[6], 3))
#                 statsData['70th Percentile'].append(round(percentiles[7], 3))
#                 statsData['75th Percentile'].append(round(percentiles[8], 3))
#                 statsData['80th Percentile'].append(round(percentiles[9], 3))
#                 statsData['90th Percentile'].append(round(percentiles[10], 3))
#                 statsData['Max'].append(round(stats['max'], 3))
#                 statsData['1 Std Dev Above'].append(round(mean + std_dev, 3))
#                 statsData['1 Std Dev Below'].append(round(mean - std_dev, 3))
#                 statsData['2 Std Dev Above'].append(round(mean + 2 * std_dev, 3))
#                 statsData['2 Std Dev Below'].append(round(mean - 2 * std_dev, 3))

#             return pd.DataFrame(statsData)

#         # Main execution
#         combinedStats = combineChunkStats(input.chunkStatsData)
#         familyStatisticsDf = generateFamilyStatistics(combinedStats)

#         # Save family statistics to CSV
#         familyStatisticsDf.to_csv(output.familyStats, index=False)

###

###
# Version 3.2

# rule eval_apply_reverse_seqs:
#     input:
#         data=join(
#             "output",
#             "vector"
#             if config["learnapp"]["fragmentation"] == False
#             else "vector_frag",
#             "{nb}.npz",
#         ),
#         annotation=expand("{an}", an=annotFiles),
#         compareAssociations=join("output", "learn", "kmer-counts-total.csv"),
#     output:
#         chunkStats=join(
#             "output",
#             "eval_apply_reversed"
#             if config["learnapp"]["fragmentation"] == False
#             else "eval_apply_frag",
#             "chunk-stats-{nb}.csv",
#         ),
#     log:
#         join(outDir, "eval_apply_reversed", "log", "{nb}.log"),
#     run:
#         start_time = datetime.now()
#         with open(log[0], "a") as f:
#             f.write(f"start time:\t{start_time}\n")

#         class KmerCompare:
#             """
#             Initializes the KmerCompare object.

#             This object is designed to compare kmer counts with provided annotations.

#             Attributes:
#                 compareAssociations (str): Path to a CSV file containing kmer counts totals matrix.
#                 annotationFiles (list): List of paths to files containing sequence annotations.
#                 inputData (str): Path to input data for kmer analysis.
#                 outputPath (str): Path to save the result.
#                 seqAnnot (dict): Dictionary mapping sequence IDs to annotations.
#                 kmerCountTotals (DataFrame): DataFrame of kmer count totals from compareAssociations.
#                 kmerList (list): List of unique kmers found in input data.
#                 bin_edges (ndarray): Array of bin edges for histograms.
#                 num_bins (int): Number of bins for histograms.
#             """

#             def __init__(
#                 self, compareAssociations, annotationFiles, inputData, outputPath
#             ):
#                 self.compareAssociations = compareAssociations
#                 self.annotationFiles = annotationFiles
#                 self.inputData = inputData
#                 self.outputPath = outputPath
#                 self.seqAnnot = {}
#                 self.kmerCountTotals = None
#                 self.kmerList = []
#                 self.num_bins = 100
#                 self.bin_edges = np.linspace(-1, 1, self.num_bins + 1)  # Bins from -1 to 1

#             def generateInputs(self):
#                 """
#                 Generates the necessary inputs for comparison.

#                 Loads kmer counts and annotations into appropriate data structures.
#                 """
#                 self.kmerCountTotals = pd.read_csv(
#                     str(self.compareAssociations),
#                     index_col="__index_level_0__",
#                     header=0,
#                     engine="c",
#                 )
#                 for f in self.annotationFiles:
#                     annotation_df = pd.read_table(f)
#                     seqs = annotation_df["id"].tolist()
#                     anns = annotation_df["Family"].tolist()
#                     for seqid, ann in zip(seqs, anns):
#                         self.seqAnnot[seqid] = ann

#             def executeAll(self, config):
#                 """
#                 Executes all the comparison steps in sequence.

#                 This includes:
#                     1. Loading inputs.
#                     2. Processing each sequence in the input data.
#                     3. Generating kmer counts for reversed sequences.
#                     4. Calculating cosine similarity and updating family statistics.
#                     5. Writing per-chunk family statistics.
#                 """
#                 self.generateInputs()
#                 # Initialize chunkStats
#                 chunkStats = {}
#                 # Load the input data
#                 kmerList, df = skm.io.load_npz(self.inputData)
#                 self.kmerList = kmerList[0]
#                 seqids = df["sequence_id"]
#                 sequences = df["sequence"]
#                 k_len = len(self.kmerList[0])

#                 for seq_id, seq in zip(seqids, sequences):
#                     v = seq[::-1]
#                     # Generate kmer counts
#                     kCounts = {}
#                     items = [v[i : i + k_len] for i in range(len(v) - k_len + 1)]
#                     for kmer in items:
#                         kCounts[kmer] = kCounts.get(kmer, 0) + 1
#                     # Create kmerCounts dataframe
#                     kmerCounts = pd.DataFrame([kCounts], index=[seq_id])
#                     kmerCounts = kmerCounts.reindex(columns=self.kmerList, fill_value=0)
#                     # Match format
#                     kmerCounts = self.matchkmerCountsFormat(kmerCounts)
#                     # Compute cosine similarities
#                     cosine_scores = self.calculateCosineSimilarity(kmerCounts)
#                     # Update chunkStats
#                     for family, score in cosine_scores.items():
#                         if family not in chunkStats:
#                             chunkStats[family] = {
#                                 'count': 0,
#                                 'sum': 0.0,
#                                 'sumSqr': 0.0,
#                                 'min': np.inf,
#                                 'max': -np.inf,
#                                 'hist_counts': np.zeros(self.num_bins, dtype=int)
#                             }
#                         stats = chunkStats[family]
#                         stats['count'] += 1
#                         stats['sum'] += score
#                         stats['sumSqr'] += score ** 2
#                         stats['min'] = min(stats['min'], score)
#                         stats['max'] = max(stats['max'], score)
#                         # Update histogram
#                         bin_index = np.searchsorted(self.bin_edges, score, side='right') - 1
#                         if 0 <= bin_index < self.num_bins:
#                             stats['hist_counts'][bin_index] += 1

#                 # Write chunkStats to output file
#                 self.writeChunkStats(chunkStats)

#             def matchkmerCountsFormat(self, kmerCounts):
#                 """
#                 Matches the format of the provided kmer counts DataFrame to the format of the
#                 comparison data. Ensures columns align correctly.
#                 """
#                 kmerCounts = kmerCounts.reindex(
#                     columns=self.kmerCountTotals.columns, fill_value=0
#                 )
#                 return kmerCounts

#             def calculateCosineSimilarity(self, kmerCounts):
#                 """
#                 Calculates the cosine similarity between the input kmer counts and comparison data.

#                 Returns:
#                     Series: A Series of cosine similarity scores.
#                 """
#                 cosine_scores = sklearn.metrics.pairwise.cosine_similarity(
#                     self.kmerCountTotals, kmerCounts
#                 ).flatten()
#                 cosine_scores_series = pd.Series(
#                     cosine_scores, index=self.kmerCountTotals.index
#                 )
#                 return cosine_scores_series

#             def writeChunkStats(self, chunkStats):
#                 """
#                 Writes the per-chunk family statistics to a CSV file.
#                 """
#                 data = []
#                 for family, stats in chunkStats.items():
#                     row = {
#                         'Family': family,
#                         'count': stats['count'],
#                         'sum': stats['sum'],
#                         'sumSqr': stats['sumSqr'],
#                         'min': stats['min'],
#                         'max': stats['max'],
#                     }
#                     # Add histogram counts
#                     for i in range(self.num_bins):
#                         row[f'hist_bin_{i}'] = stats['hist_counts'][i]
#                     data.append(row)
#                 df = pd.DataFrame(data)
#                 df.to_csv(self.outputPath, index=False)

#         analysis = KmerCompare(
#             input.compareAssociations, input.annotation, input.data, output.chunkStats
#         )
#         analysis.executeAll(config)
#         skm.utils.log_runtime(log[0], start_time)

# rule reverseDecoy_evaluations:
#     input:
#         chunkStatsData=expand(
#             join(
#                 "output",
#                 "eval_apply_reversed"
#                 if config["learnapp"]["fragmentation"] == False
#                 else "eval_apply_frag",
#                 "chunk-stats-{nb}.csv",
#             ),
#             nb=FAS,
#         ),
#     output:
#         familyStats="output/eval_conf/family_summary_stats.csv",
#     run:
#         import pandas as pd
#         import numpy as np

#         def combineChunkStats(chunkStatsFiles, num_bins):
#             combinedStats = {}
#             for filename in chunkStatsFiles:
#                 df = pd.read_csv(filename)
#                 for index, row in df.iterrows():
#                     family = row['Family']
#                     if family not in combinedStats:
#                         combinedStats[family] = {
#                             'count': 0,
#                             'sum': 0.0,
#                             'sumSqr': 0.0,
#                             'min': np.inf,
#                             'max': -np.inf,
#                             'hist_counts': np.zeros(num_bins, dtype=int)
#                         }
#                     stats = combinedStats[family]
#                     stats['count'] += row['count']
#                     stats['sum'] += row['sum']
#                     stats['sumSqr'] += row['sumSqr']
#                     stats['min'] = min(stats['min'], row['min'])
#                     stats['max'] = max(stats['max'], row['max'])
#                     # Update histogram counts
#                     for i in range(num_bins):
#                         stats['hist_counts'][i] += row[f'hist_bin_{i}']
#             return combinedStats

#         def generateFamilyStatistics(combinedStats, bin_edges, num_bins):
#             """
#             Generates statistics for each family using the combined statistics.

#             Args:
#                 combinedStats (dict): Combined statistics for each family.

#             Returns:
#                 pd.DataFrame: DataFrame containing the statistics for each family.
#             """
#             statsData = {
#                 'Family': [],
#                 'Mean': [],
#                 'Std Dev': [],
#                 'Min': [],
#                 '10th Percentile': [],
#                 '20th Percentile': [],
#                 '25th Percentile': [],
#                 '30th Percentile': [],
#                 '40th Percentile': [],
#                 'Median': [],
#                 '60th Percentile': [],
#                 '70th Percentile': [],
#                 '75th Percentile': [],
#                 '80th Percentile': [],
#                 '90th Percentile': [],
#                 'Max': [],
#                 '1 Std Dev Above': [],
#                 '1 Std Dev Below': [],
#                 '2 Std Dev Above': [],
#                 '2 Std Dev Below': [],
#             }

#             for family, stats in combinedStats.items():
#                 n = stats['count']
#                 sum_ = stats['sum']
#                 sumSqr = stats['sumSqr']
#                 mean = sum_ / n
#                 variance = (
#                     (sumSqr - (sum_ ** 2) / n) / (n - 1) if n > 1 else 0.0
#                 )
#                 std_dev = np.sqrt(variance)

#                 # Approximate percentiles from histogram
#                 cumsum = np.cumsum(stats['hist_counts'])
#                 total = cumsum[-1]
#                 percentiles_values = {}
#                 percentiles = [10, 20, 25, 30, 40, 50, 60, 70, 75, 80, 90]
#                 for p in percentiles:
#                     idx = np.searchsorted(cumsum, total * p / 100)
#                     if idx >= num_bins:
#                         idx = num_bins - 1
#                     percentiles_values[p] = (bin_edges[idx] + bin_edges[idx + 1]) / 2

#                 statsData['Family'].append(family)
#                 statsData['Mean'].append(round(mean, 3))
#                 statsData['Std Dev'].append(round(std_dev, 3))
#                 statsData['Min'].append(round(stats['min'], 3))
#                 statsData['10th Percentile'].append(round(percentiles_values[10], 3))
#                 statsData['20th Percentile'].append(round(percentiles_values[20], 3))
#                 statsData['25th Percentile'].append(round(percentiles_values[25], 3))
#                 statsData['30th Percentile'].append(round(percentiles_values[30], 3))
#                 statsData['40th Percentile'].append(round(percentiles_values[40], 3))
#                 statsData['Median'].append(round(percentiles_values[50], 3))
#                 statsData['60th Percentile'].append(round(percentiles_values[60], 3))
#                 statsData['70th Percentile'].append(round(percentiles_values[70], 3))
#                 statsData['75th Percentile'].append(round(percentiles_values[75], 3))
#                 statsData['80th Percentile'].append(round(percentiles_values[80], 3))
#                 statsData['90th Percentile'].append(round(percentiles_values[90], 3))
#                 statsData['Max'].append(round(stats['max'], 3))
#                 statsData['1 Std Dev Above'].append(round(mean + std_dev, 3))
#                 statsData['1 Std Dev Below'].append(round(mean - std_dev, 3))
#                 statsData['2 Std Dev Above'].append(round(mean + 2 * std_dev, 3))
#                 statsData['2 Std Dev Below'].append(round(mean - 2 * std_dev, 3))

#             return pd.DataFrame(statsData)

#         # Main execution
#         num_bins = 100
#         bin_edges = np.linspace(-1, 1, num_bins + 1)
#         combinedStats = combineChunkStats(input.chunkStatsData, num_bins)
#         familyStatisticsDf = generateFamilyStatistics(combinedStats, bin_edges, num_bins)

#         # Save family statistics to CSV
#         familyStatisticsDf.to_csv(output.familyStats, index=False)


###
###
###


rule eval_apply_sequences:
    input:
        data=join(
            "output",
            "vector"
            if config["learnapp"]["fragmentation"] == False
            else "vector_frag",
            "{nb}.npz",
        ),
        annotation=expand("{an}", an=annotFiles),
        compareAssociations=join("output", "learn", "kmer-counts-total.csv"),
    output:
        apply=join(
            "output",
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
                            (x + "_unknown_" + str(count))
                        ] = self.seqKmerdict.pop(seqid)
                    else:
                        self.seqKmerdict[
                            (self.seqAnnot[x] + "_known_" + str(count))
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
                top_1_indices = np.argsort(-finalMatrixWithScores.values, axis=1)[:, :1]  # Changed to top 1
                mask = np.zeros_like(finalMatrixWithScores.values, dtype=bool)

                for i, indexOne in enumerate(top_1_indices):  # Only one index now
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
                "output",
                "eval_apply_sequences"
                if config["learnapp"]["fragmentation"] == False
                else "eval_apply_frag",
                "seq-annotation-scores-{nb}.csv.gz",
            ),
            nb=FAS,
        ),
        baseConfidence=expand("{bc}", bc=baseConfidence),
        reverseDecoyStats="output/eval_conf/family_summary_stats.csv",
    output:
        evalGlob="output/eval_conf/global-confidence-scores.csv",
    params:
        modifier=config["learnapp"]["conf_weight_modifier"],
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
                threshold_dict = dict(zip(thresholds_df.Family, thresholds_df[threshold_type]))
                self.threshold_dict = threshold_dict  # Store for later use

                # Apply the selection method
                predictions, deltas, topTwo = self.apply_selectionMethod(seqAnnScores)
                
                result = seqAnnScores.index.tolist()  # Assuming the actual labels are in the index

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

            def apply_selectionMethod(self, seqAnnScores):
                selectionMethod = config["learnapp"]["selection"]
                threshold_type = config["learnapp"]["threshold"]

                if selectionMethod == 'top_hit_no_threshold':

                    #Doesn't work for PFam for some reason.
                    # Method 0: Select top hit without threshold
                    # Get the top two scores and indices
                    # print(f"seqAnnScores.values: {seqAnnScores.values}")
                    # # Use numpy to partition the values, get the two largest without fully sorting
                    # partitioned_vals = np.partition(seqAnnScores.values, -2, axis=1)
                    # print(f"partitioned_vals: {partitioned_vals}")

                    # # Get the indices of the top 2 values
                    # topTwoIdx = np.argpartition(seqAnnScores.values, -2, axis=1)[:, -2:]
                    # print(f"topTwoIdx: {topTwoIdx}")

                    # # Sort the indices to get the max and second max in the correct order
                    # sortedTopTwoIdx = np.argsort(seqAnnScores.values[np.arange(seqAnnScores.shape[0])[:, None], topTwoIdx], axis=1)[:, ::-1]
                    # print(f"sortedTopTwoIdx: {sortedTopTwoIdx}")

                    # # Get the actual top 2 values
                    # top_2_vals = np.take_along_axis(seqAnnScores.values, topTwoIdx[np.arange(seqAnnScores.shape[0])[:, None], sortedTopTwoIdx], axis=1)
                    # print(f"top_2_vals: {top_2_vals}")

                    # # Get the corresponding column headers for the top 2 values
                    # topTwoHeaders = np.array(seqAnnScores.columns)[topTwoIdx[np.arange(seqAnnScores.shape[0])[:, None], sortedTopTwoIdx]]
                    # print(f"topTwoHeaders: {topTwoHeaders}")

                    # # Create a new DataFrame with the two highest values and their respective column headers
                    # keyValsDf = pd.DataFrame({
                    #     'keyValueOne': top_2_vals[:, 0],  # Largest value in each row
                    #     'keyValueOneHeader': topTwoHeaders[:, 0],  # Column name of the largest value
                    #     'keyValueTwo': top_2_vals[:, 1],  # Second largest value in each row
                    #     'keyValueTwoHeader': topTwoHeaders[:, 1]   # Column name of the second largest value
                    # })
                    

                    # More Memory Intensive
                    def get_topTwo(row):
                        top2 = row.nlargest(2)
                        return pd.Series({
                            'keyValueOne': top2.iloc[0] if len(top2) > 0 else np.nan,
                            'keyValueOneHeader': top2.index[0] if len(top2) > 0 else np.nan,
                            'keyValueTwo': top2.iloc[1] if len(top2) > 1 else np.nan,
                            'keyValueTwoHeader': top2.index[1] if len(top2) > 1 else np.nan
                        })

                    # Apply the function to each row of your DataFrame
                    keyValsDf = seqAnnScores.apply(get_topTwo, axis=1)

                    # Display the resulting DataFrame
                    print(keyValsDf)

                    predictions = keyValsDf['keyValueOneHeader'].tolist()  # List of largest value headers
                    deltas = (keyValsDf['keyValueOne'] - keyValsDf['keyValueTwo']).tolist()  # Difference between largest and second largest values
                    topTwo = keyValsDf[['keyValueOne', 'keyValueTwo']]  # DataFrame with keyValueOne and keyValueTwo


                    #Let see... also memory fail
                    # seqAnnScores_filled = seqAnnScores.fillna(-np.inf)

                    # # Get the numpy array and column names
                    # values = seqAnnScores_filled.values
                    # column_names = seqAnnScores_filled.columns.to_numpy()

                    # # Get the indices that would sort the values in descending order along each row
                    # sorted_indices = np.argsort(-values, axis=1)

                    # # Prepare row indices
                    # row_indices = np.arange(values.shape[0])[:, None]

                    # # Get the top two indices for each row
                    # top_n = 2  # Number of top values to retrieve
                    # top_indices = sorted_indices[:, :top_n]

                    # # Get the top two values and their corresponding headers
                    # top_values = values[row_indices, top_indices]
                    # top_headers = column_names[top_indices]

                    # # Convert -infinity back to NaN for clarity
                    # top_values[top_values == -np.inf] = np.nan

                    # # Create a DataFrame for the top two values and headers
                    # keyValsDf = pd.DataFrame({
                    #     'keyValueOne': top_values[:, 0],
                    #     'keyValueOneHeader': top_headers[:, 0],
                    #     'keyValueTwo': top_values[:, 1],
                    #     'keyValueTwoHeader': top_headers[:, 1]
                    # })

                    # # Display the resulting DataFrame
                    # print(keyValsDf)

                    # # Extract predictions, deltas, and topTwo as per your original code
                    # predictions = keyValsDf['keyValueOneHeader'].tolist()
                    # deltas = (keyValsDf['keyValueOne'] - keyValsDf['keyValueTwo']).tolist()
                    # topTwo = keyValsDf[['keyValueOne', 'keyValueTwo']]
                    return predictions, deltas, topTwo

                elif selectionMethod == 'top_hit_with_threshold':
                    # # Method 1: Top hit above threshold
                    # # Map thresholds to the columns (families)
                    # thresholds = seqAnnScores.columns.to_series().map(self.threshold_dict)

                    # filteredScores = seqAnnScores.where(seqAnnScores >= thresholds, np.nan)

                    # # Method: Select top two scores above threshold
                    # # Use numpy to partition the values, get the two largest without fully sorting (ignoring NaNs)
                    # partitioned_vals = np.partition(np.nan_to_num(filteredScores.values, nan=-np.inf), -2, axis=1)

                    # # Get the indices of the top 2 values
                    # topTwoIdx = np.argpartition(np.nan_to_num(filteredScores.values, nan=-np.inf), -2, axis=1)[:, -2:]

                    # # Sort the indices to get the max and second max in the correct order
                    # sortedTopTwoIdx = np.argsort(filteredScores.values[np.arange(filteredScores.shape[0])[:, None], topTwoIdx], axis=1)[:, ::-1]

                    # # Get the actual top 2 values
                    # top_2_vals = np.take_along_axis(filteredScores.values, topTwoIdx[np.arange(filteredScores.shape[0])[:, None], sortedTopTwoIdx], axis=1)

                    # # Get the corresponding column headers for the top 2 values
                    # topTwoHeaders = np.array(filteredScores.columns)[topTwoIdx[np.arange(filteredScores.shape[0])[:, None], sortedTopTwoIdx]]

                    # # Create a new DataFrame with the two highest values and their respective column headers
                    # keyValsDf = pd.DataFrame({
                    #     'keyValueOne': top_2_vals[:, 0],  # Largest value in each row
                    #     'keyValueOneHeader': topTwoHeaders[:, 0],  # Column name of the largest value
                    #     'keyValueTwo': top_2_vals[:, 1],  # Second largest value in each row
                    #     'keyValueTwoHeader': topTwoHeaders[:, 1]   # Column name of the second largest value
                    # })
                    # print(f"keyValsDf: {keyValsDf}")
                    # predictions = keyValsDf['keyValueOneHeader'].tolist()  # List of largest value headers
                    # deltas = (keyValsDf['keyValueOne'] - keyValsDf['keyValueTwo']).tolist()  # Difference between largest and second largest values
                    # topTwo = keyValsDf[['keyValueOne', 'keyValueTwo']]  # DataFrame with keyValueOne and keyValueTwo

                    # return predictions, deltas, topTwo
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
                    predictions = keyValsDf['keyValueOneHeader'].tolist()  # List of largest value headers
                    deltas = (keyValsDf['keyValueOne'] - keyValsDf['keyValueTwo']).tolist()  # Difference between largest and second largest values
                    topTwo = keyValsDf[['keyValueOne', 'keyValueTwo']]  # DataFrame with keyValueOne and keyValueTwo

                    return predictions, deltas, topTwo


                # Delta compares value from threshold  - Lets call this Method: Val - Thresh
                elif selectionMethod == 'greatest_distance_from_threshold_dt':
                    thresholds = seqAnnScores.columns.to_series().map(self.threshold_dict)
                    filteredScores = seqAnnScores.where(seqAnnScores >= thresholds, np.nan)

                    # Method: Select top score above threshold
                    # Get the top value and its corresponding column
                    top_vals = filteredScores.max(axis=1)
                    top_headers = filteredScores.idxmax(axis=1)
                    print(f"thresholds: {thresholds}")
                    print(f"top_vals: {top_vals}")
                    print(f"top_headers: {top_headers}")
                    # Create a new DataFrame with the top value, threshold, and their respective column headers
                    keyValsDf = pd.DataFrame({
                        'keyValueOne': top_vals,  # Largest value in each row
                        'keyValueOneHeader': top_headers,  # Column name of the largest value
                        'keyValueTwo': top_headers.map(thresholds).values,  # Threshold value for the largest value
                        'keyValueTwoHeader': top_headers  # Column name of the threshold value
                    })

                    predictions = keyValsDf['keyValueOneHeader'].tolist()  # List of largest value headers
                    deltas = (keyValsDf['keyValueOne'] - keyValsDf['keyValueTwo']).tolist()  # Difference between largest value and threshold
                    topTwo = keyValsDf[['keyValueOne', 'keyValueTwo']] 

                    return predictions, deltas, topTwo


                # # # Maybe more accurate.  Delta compares top greatest distances
                elif selectionMethod == 'greatest_distance_from_threshold_tt':
                    thresholds = seqAnnScores.columns.to_series().map(self.threshold_dict)
                    
                    # Compute the distances from thresholds
                    distances = seqAnnScores.subtract(thresholds, axis=1)
                    
                    # Filter out distances less than 0 (scores below threshold)
                    filtered_distances = distances.where(distances >= 0, np.nan)
                    
                    # Method: Select top two distances
                    # Use numpy to partition the values, get the two largest distances without fully sorting (ignoring NaNs)
                    partitioned_distances = np.partition(
                        np.nan_to_num(filtered_distances.values, nan=-np.inf), -2, axis=1
                    )
                    
                    # Get the indices of the top 2 distances
                    topTwoIdx = np.argpartition(
                        np.nan_to_num(filtered_distances.values, nan=-np.inf), -2, axis=1
                    )[:, -2:]
                    
                    # Sort the indices to get the max and second max in the correct order
                    sortedTopTwoIdx = np.argsort(
                        filtered_distances.values[
                            np.arange(filtered_distances.shape[0])[:, None], topTwoIdx
                        ],
                        axis=1,
                    )[:, ::-1]
                    
                    # Get the actual top 2 distances
                    topTwoDistances = np.take_along_axis(
                        filtered_distances.values,
                        topTwoIdx[np.arange(filtered_distances.shape[0])[:, None], sortedTopTwoIdx],
                        axis=1,
                    )
                    
                    # Get the corresponding scores
                    topTwoScores = np.take_along_axis(
                        seqAnnScores.values,
                        topTwoIdx[np.arange(filtered_distances.shape[0])[:, None], sortedTopTwoIdx],
                        axis=1,
                    )
                    
                    # Get the corresponding column headers for the top 2 distances
                    topTwoHeaders = np.array(filtered_distances.columns)[
                        topTwoIdx[np.arange(filtered_distances.shape[0])[:, None], sortedTopTwoIdx]
                    ]
                    
                    # Get thresholds for the top 2 headers
                    # Flatten the headers to 1D array
                    flattened_headers = topTwoHeaders.flatten()
                    
                    # Retrieve thresholds and reshape back to (n_samples, 2)
                    flattened_thresholds = thresholds[flattened_headers].values
                    topTwoThresholds = flattened_thresholds.reshape(topTwoHeaders.shape)
                    
                    # Create a new DataFrame with the top values, distances, thresholds, and their respective column headers
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
                    
                    # Calculate deltas as the difference between the distances from thresholds
                    deltas = (keyValsDf['keyValueOneDistance'] - keyValsDf['keyValueTwoDistance']).tolist()
                    
                    predictions = keyValsDf['keyValueOneHeader'].tolist()
                    topTwo = keyValsDf[['keyValueOne', 'keyValueTwo']]
                    
                    return predictions, deltas, topTwo

        
                elif selectionMethod == 'combined_method_with_dynamic_delta':
                    # Method 4: Combined method with dynamic delta calculation
                    # Set weights
                    weight_top = config["learnapp"].get("weight_top", 0.5)
                    weight_distance = config["learnapp"].get("weight_distance", 0.5)

                    # Map thresholds to the columns (families)
                    thresholds = seqAnnScores.columns.to_series().map(self.threshold_dict)

                    # Calculate distances
                    distances = seqAnnScores - thresholds

                    # Only consider positive distances
                    positiveDistances = distances.where(distances >= 0, np.nan)

                    # Calculate combined scores
                    combinedScores = (seqAnnScores * weight_top) + (distances * weight_distance)

                    # Only consider combined scores where distances are positive
                    combinedScores = combinedScores.where(positiveDistances.notna(), np.nan)

                    # Select the top family based on combined scores
                    topCombinedScores = combinedScores.max(axis=1)
                    predictions = combinedScores.idxmax(axis=1)
                    predictions = predictions.where(~topCombinedScores.isna(), None)

                    # Prepare data structures for delta calculation
                    deltas = []
                    topTwoList = []

                    # Precompute top families and scores from methods 1/2 and method 3

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
                            # For topTwo, use top two scores from method 1/2
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

                ###first attempt 90%
                elif selectionMethod == 'combined_topTwo_scores':
                    # Method 4: Combine methods 2 and 3
                    # Set weights
                    weight_top = config["learnapp"].get("weight_top", 0.5)
                    weight_distance = config["learnapp"].get("weight_distance", 0.5)

                    # Map thresholds to the columns (families)
                    thresholds = seqAnnScores.columns.to_series().map(self.threshold_dict)

                    # Calculate distances
                    distances = seqAnnScores - thresholds

                    # Only consider positive distances
                    positiveDistances = distances.where(distances >= 0, np.nan)

                    # Calculate combined scores
                    combinedScores = (seqAnnScores * weight_top) + (distances * weight_distance)

                    # Only consider combined scores where distances are positive
                    combinedScores = combinedScores.where(positiveDistances.notna(), np.nan)

                    # Method: Select top two combined scores
                    # Use numpy to partition the values, get the two largest combined scores without fully sorting (ignoring NaNs)
                    partitioned_combinedScores = np.partition(
                        np.nan_to_num(combinedScores.values, nan=-np.inf), -2, axis=1
                    )

                    # Get the indices of the top 2 combined scores
                    topTwoIdx = np.argpartition(
                        np.nan_to_num(combinedScores.values, nan=-np.inf), -2, axis=1
                    )[:, -2:]

                    # Sort the indices to get the max and second max in the correct order
                    sortedTopTwoIdx = np.argsort(
                        combinedScores.values[
                            np.arange(combinedScores.shape[0])[:, None], topTwoIdx
                        ],
                        axis=1
                    )[:, ::-1]

                    # Get the actual top 2 combined scores
                    topTwoCombinedScores = np.take_along_axis(
                        combinedScores.values,
                        topTwoIdx[np.arange(combinedScores.shape[0])[:, None], sortedTopTwoIdx],
                        axis=1,
                    )

                    # Get the corresponding original scores
                    topTwoScores = np.take_along_axis(
                        seqAnnScores.values,
                        topTwoIdx[np.arange(combinedScores.shape[0])[:, None], sortedTopTwoIdx],
                        axis=1,
                    )

                    # Get the corresponding distances
                    topTwoDistances = np.take_along_axis(
                        distances.values,
                        topTwoIdx[np.arange(combinedScores.shape[0])[:, None], sortedTopTwoIdx],
                        axis=1,
                    )

                    # Get the corresponding thresholds
                    topTwoThresholds = np.take_along_axis(
                        thresholds.values[np.newaxis, :],
                        topTwoIdx[np.arange(combinedScores.shape[0])[:, None], sortedTopTwoIdx],
                        axis=1,
                    )

                    # Get the corresponding column headers for the top 2 combined scores
                    topTwoHeaders = np.array(combinedScores.columns)[
                        topTwoIdx[np.arange(combinedScores.shape[0])[:, None], sortedTopTwoIdx]
                    ]

                    # Create a new DataFrame with the top values, distances, thresholds, and their respective column headers
                    keyValsDf = pd.DataFrame({
                        'keyValueOne': topTwoCombinedScores[:, 0],
                        'keyValueOneHeader': topTwoHeaders[:, 0],
                        'keyValueTwo': topTwoCombinedScores[:, 1],
                        'keyValueTwoHeader': topTwoHeaders[:, 1],
                        'keyValueOneScore': topTwoScores[:, 0],
                        'keyValueOneDistance': topTwoDistances[:, 0],
                        'keyValueOneThreshold': topTwoThresholds[:, 0],
                        'keyValueTwoScore': topTwoScores[:, 1],
                        'keyValueTwoDistance': topTwoDistances[:, 1],
                        'keyValueTwoThreshold': topTwoThresholds[:, 1],
                    })

                    # Predictions are the top headers
                    predictions = keyValsDf['keyValueOneHeader'].tolist()

                    # Deltas are the difference between top two combined scores
                    deltas = (keyValsDf['keyValueOne'] - keyValsDf['keyValueTwo']).tolist()

                    # topTwo DataFrame with keyValueOne and keyValueTwo
                    topTwo = keyValsDf[['keyValueOne', 'keyValueTwo']]

                    return predictions, deltas, topTwo
                ###

                # elif selectionMethod == 'combined_score_with_threshold':
                #     # Method 3: Balanced distance
                #     # Set weights
                #     weight_top = config["learnapp"].get("weight_top", 0.5)
                #     weight_distance = config["learnapp"].get("weight_distance", 0.5)

                #     # Map thresholds to the columns (families)
                #     thresholds = seqAnnScores.columns.to_series().map(self.threshold_dict)
                #     # Calculate distances
                #     distances = seqAnnScores - thresholds
                #     # Only consider positive distances
                #     positiveDistances = distances.where(distances > 0, np.nan)
                #     # Calculate combined scores
                #     combinedScores = (seqAnnScores * weight_top) + (distances * weight_distance)
                #     # Only consider combined scores where distances are positive
                #     combinedScores = combinedScores.where(positiveDistances.notna(), np.nan)
                #     # Select the family with the highest combined score
                #     topCombinedScores = combinedScores.max(axis=1)
                #     predictions = combinedScores.idxmax(axis=1)
                #     predictions = predictions.where(~topCombinedScores.isna(), None)

                #     # Delta is the combined score
                #     deltas = topCombinedScores.round(2)

                #     return predictions, deltas

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
                

                # print("diffDataframe:\n", diffDataframe.head())
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
                Computes the ratio distribution for True and False total distributions.

                Args:
                    trueTotalDist (Series): Total distribution for True predictions.
                    false_total_dist (Series): Total distribution for False predictions.

                Returns:
                    Series: Computed ratio distribution.
                """
                # Calculate the ratio
                ratioTotalDist = trueTotalDist / (trueTotalDist + false_total_dist)
                
                # Create a new index from 0 to 1 with .01 increments
                newIndex = pd.Index([round(i, 2) for i in pd.np.arange(0, 1.01, 0.01)], name="ratioTotalDist")
                ratioTotalDist = ratioTotalDist.reindex(newIndex)
                
                # Interpolate the missing values linearly and fill nan values with 0.
                ratioTotalDist = ratioTotalDist.interpolate(method="linear")
                ratioTotalDist.fillna(0, inplace=True)
                return ratioTotalDist

            def checkConfidenceMerge(self, newRatioDist):
                """
                Merge with base confidence file if available.

                Args:
                    newRatioDist (Series): Total distribution for T/(T+F) predictions.

                Returns:
                    DataFrame: Updated ratio distribution.
                """
                sumSeries = (
                    self.trueRunningCrosstab.sum() + self.falseRunningCrosstab.sum()
                )
                currentWeight = (
                    self.trueRunningCrosstab.values.sum()
                    + self.falseRunningCrosstab.values.sum()
                )

                updatedData = pd.DataFrame(newRatioDist, columns=["confidence"])
                updatedData = updatedData.assign(
                    weight=[currentWeight] * len(updatedData), sum=sumSeries
                )
                updatedData.index.name = 'Difference'

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
                    priorWeightedScore = (
                        priorConf["confidence"] * priorConf["weight"]
                    )
                    currentWeighted_score = (
                        updatedData["confidence"] * weightedCurrent
                    )

                    updatedConfidence = (
                        priorWeightedScore + currentWeighted_score
                    ) / totalWeight

                    updatedData = pd.DataFrame({"confidence": updatedConfidence})

                    updatedData = updatedData.assign(
                        weight=outWeight,
                        sum=sumSeries + priorConf["sum"],
                    )

                    print(f"Final Confidence Data\n{updatedData}")
                else:
                    print(
                        "Base confidence file not found or multiple files present. Only one file is allowed in baseConfidence."
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
                ratioDist = self.computeRatioDistribution(trueDist, falseDist)
                ratioDist = self.checkConfidenceMerge(ratioDist)
                globalConfidence = self.checkConfidenceMerge(ratioDist)
                # globalConfidence = self.smooth_confidence(globalConfidence)
                self.saveResults(
                    globalConfidence,
                    self.outputGlob
                )

        evaluator = Evaluator(
            input.evalApplyData,
            output.evalGlob,
            input.reverseDecoyStats,
            params.modifier,  # Pass the modifier to the Evaluator
            input.baseConfidence,
        )
        evaluator.executeAll()

        skm.utils.log_runtime(log[0], start_time)
