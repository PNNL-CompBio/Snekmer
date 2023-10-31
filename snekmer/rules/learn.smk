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

# Note:
# Pyarrow installed via "conda install -c conda-forge pyarrow"
# collect all fasta-like files, unzipped filenames, and basenames
input_dir = (
    "input"
    if (("input_dir" not in config) or (str(config["input_dir"]) == "None"))
    else config["input_dir"]
)
input_files = glob(join(input_dir, "*"))
# base_file = glob(join(input_dir,"base" "*"))
zipped = [fa for fa in input_files if fa.endswith(".gz")]
unzipped = [
    fa.rstrip(".gz")
    for fa, ext in product(input_files, config["input_file_exts"])
    if fa.rstrip(".gz").endswith(f".{ext}")
]
annot_files = glob(join("annotations", "*.ann"))
base_counts = glob(join("base", "counts", "*.csv"))
base_confidence = glob(join("base", "confidence", "*.csv"))
# map extensions to basename (basename.ext.gz -> {basename: ext})
UZ_MAP = {
    skm.utils.split_file_ext(f)[0]: skm.utils.split_file_ext(f)[1] for f in zipped
}
FA_MAP = {
    skm.utils.split_file_ext(f)[0]: skm.utils.split_file_ext(f)[1] for f in unzipped
}

# get unzipped filenames
UZS = [f"{f}.{ext}" for f, ext in UZ_MAP.items()]
# isolate basenames for all files
FAS = list(FA_MAP.keys())
# parse any background files
bg_files = glob(join(input_dir, "background", "*"))
if len(bg_files) > 0:
    bg_files = [skm.utils.split_file_ext(basename(f))[0] for f in bg_files]
NON_BGS, BGS = [f for f in FAS if f not in bg_files], bg_files
# terminate with error if invalid alphabet specified
skm.alphabet.check_valid(config["alphabet"])
# define output directory (helpful for multiple runs)
out_dir = skm.io.define_output_dir(
    config["alphabet"], config["k"], nested=config["nested_output"]
)


output_prefixes = (
    ["vector"] if not config["learnapp"]["fragmentation"] else ["vector", "vector_frag"]
)


# define output files to be created by snekmer
rule all:
    input:
        expand(join("output", "vector", "{nb}.npz"), nb=FAS),
        expand(join("output", "learn", "kmer-counts-{nb}.csv"), nb=FAS),
        join("output", "learn", "kmer-counts-total.csv"),
        expand(join("output", "fragmented", "{nb}.fasta"), nb=FAS)
        if config["learnapp"]["fragmentation"]
        else [],
        expand(join("output", "vector_frag", "{nb}.npz"), nb=FAS)
        if config["learnapp"]["fragmentation"]
        else [],
        expand(
            join(
                "output",
                "eval_apply"
                if not config["learnapp"]["fragmentation"]
                else "eval_apply_frag",
                "seq-annotation-scores-{nb}.csv",
            ),
            nb=FAS,
        ),
        "output/eval_conf/confidence-matrix.csv",
        "output/eval_conf/global-confidence-scores.csv",


# if any files are gzip zipped, unzip them
use rule unzip from process with:
    output:
        unzipped=join(input_dir, "{uz}"),
        zipped=join(input_dir, "zipped", "{uz}.gz"),


if config["learnapp"]["fragmentation"]:

    rule fragmentation:
        input:
            fasta=lambda wildcards: join(
                input_dir, f"{wildcards.nb}.{FA_MAP[wildcards.nb]}"
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
                        chosen_index = random.randint(0, len(frags) - 1)
                        frags = [frags[chosen_index]]

                return frags

            with open(input.fasta, "r") as f:
                fasta_sequences = SeqIO.parse(f, "fasta")

                with open(output.fasta_out, "w") as the_file:
                    for fasta in fasta_sequences:
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
            "output" if wildcards.prefix == "vector_frag" else input_dir,
            "fragmented" if wildcards.prefix == "vector_frag" else "",
            f"{wildcards.nb}.{FA_MAP[wildcards.nb]}",
        ),
    output:
        data=join("output", "{prefix}", "{nb}.npz"),
        kmerobj=join("output", "kmerize_{prefix}", "{nb}.kmers"),
    log:
        join("output", "{prefix}_kmerize", "log", "{nb}.log"),


# WORKFLOW to learn kmer associations
# collect all seq files and generate mega-cluster
rule learn:
    input:
        data="output/vector/{nb}.npz",
        annotation=expand("{an}", an=annot_files),
    output:
        counts="output/learn/kmer-counts-{nb}.csv",
    log:
        join(out_dir, "learn", "log", "learn-{nb}.log"),
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
        seq_annot (dict): A dictionary mapping sequence IDs to their annotations.
        kmerlist (list): A list of unique kmers present in the data.
        df (DataFrame or None): DataFrame containing kmer data.
        seqids (list): A list of sequence IDs.
        kmer_totals (list): A list to store total counts of each k-mer across all sequences.
        seq_kmer_dict (dict): A dictionary mapping sequence IDs to their k-mer counts.
        annotation_counts (dict): A dictionary mapping annotations to their counts.
        total_seqs (int): The total number of sequences after filtering.
    """

            def __init__(self):
                self.annotation = []
                self.seq_annot = {}
                self.kmerlist = []
                self.df = None
                self.seqids = []
                self.kmer_totals = []
                self.seq_kmer_dict = {}
                self.annotation_counts = {}
                self.total_seqs = 0

            def load_annotations(self, input_annotation):
                """
        Load annotations from a list of provided input files.

        Args:
            input_annotation (list): List of file paths containing annotations.
        """
                for f in input_annotation:
                    self.annotation.append(pd.read_table(f))
                annotations = pd.concat(self.annotation)
                seqs = annotations["id"].tolist()
                anns = annotations["TIGRFAMs"].tolist()
                for i, seqid in enumerate(seqs):
                    self.seq_annot[seqid] = anns[i]
                self.seqs = set(seqs)

            def load_data(self, input_data):
                """
        Load and format kmer data from the provided input.

        Args:
            input_data (str): Path to the data file.
        """
                self.kmerlist, self.df = skm.io.load_npz(input_data)
                self.kmerlist = self.kmerlist[0]
                self.seqids = self.df["sequence_id"]
                for item in self.kmerlist:
                    self.kmer_totals.append(0)

            def generate_kmer_counts(self):
                """
        Generate kmer counts for sequences present in the data.
        """
                k_len = len(self.kmerlist[0])
                for i, seq in enumerate(self.seqids):
                    v = self.df["sequence"][i]
                    k_counts = self._compute_kmer_counts_for_sequence(v, k_len)
                    self.seq_kmer_dict[seq] = k_counts

            def filter_and_construct(self):
                """
        Filters sequences not present in annotations and constructs annotation counts.
        """
                self.total_seqs = len(self.seq_kmer_dict)
                for i, seqid in enumerate(list(self.seq_kmer_dict)):
                    x = re.findall(r"\|(.*?)\|", seqid)[0]
                    if x not in self.seqs:
                        del self.seq_kmer_dict[seqid]
                    else:
                        self._process_annotation_counts(seqid, x)

            def format_and_write_output(self, input_data):
                """
        Writes processed kmer counts to an output CSV file.

        Args:
            input_data (str): Path to the data file (used for naming the output file).
        """
                kmer_counts = pd.DataFrame(self.seq_kmer_dict.values())
                kmer_counts.insert(
                    0, "Annotations", self.annotation_counts.values(), True
                )

                kmer_counts_values = (
                    kmer_counts[list(kmer_counts.columns[1:])].sum(axis=1).to_list()
                )
                kmer_counts.insert(1, "Kmer Count", kmer_counts_values, True)

                self.kmer_totals[0:0] = [self.total_seqs, sum(self.kmer_totals)]
                colnames = ["Sequence count"] + ["Kmer Count"] + list(self.kmerlist)
                kmer_counts = pd.DataFrame(
                    np.insert(kmer_counts.values, 0, values=self.kmer_totals, axis=0)
                )
                kmer_counts.columns = colnames
                new_index = ["Totals"] + list(self.annotation_counts.keys())
                kmer_counts.index = new_index
                kmer_counts.replace(0, "", inplace=True)
                out_name = "output/learn/kmer-counts-" + str(input_data)[14:-4] + ".csv"
                kmer_counts.index.name = "__index_level_0__"
                kmer_counts.to_csv(out_name, index=True)
                skm.utils.log_runtime(log[0], start_time)

            def _compute_kmer_counts_for_sequence(self, v, k_len):
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
                k_counts = {}
                for j in items:
                    k_counts[j] = k_counts.get(j, 0) + 1
                store = []
                for i, item in enumerate(self.kmerlist):
                    if item in k_counts:
                        store.append(k_counts[item])
                        self.kmer_totals[i] += k_counts[item]
                    else:
                        store.append(0)
                return store

            def _process_annotation_counts(self, seqid, x):
                """
        Processes annotation counts by aggregating them based on annotation labels.

        Args:
            seqid (str): Sequence ID.
            x (str): Extracted annotation ID from seqid.
        """
                if self.seq_annot[x] not in self.seq_kmer_dict:
                    self.seq_kmer_dict[self.seq_annot[x]] = self.seq_kmer_dict.pop(
                        seqid
                    )
                else:
                    zipped_lists = zip(
                        self.seq_kmer_dict.pop(seqid),
                        self.seq_kmer_dict[self.seq_annot[x]],
                    )
                    self.seq_kmer_dict[self.seq_annot[x]] = [
                        sum(pair) for pair in zipped_lists
                    ]
                if self.seq_annot[x] not in self.annotation_counts:
                    self.annotation_counts[self.seq_annot[x]] = 1
                else:
                    self.annotation_counts[self.seq_annot[x]] += 1

            def execute_all(self, input_annotation, input_data):
                """
        Execute the entire sequence of operations in the Library process.

        Args:
            input_annotation (list): List of file paths containing annotations.
            input_data (str): Path to the data file.
        """
                self.load_annotations(input_annotation)
                self.load_data(input_data)
                self.generate_kmer_counts()
                self.filter_and_construct()
                self.format_and_write_output(input_data)


        library = Library()
        library.execute_all(input.annotation, input.data)


rule merge:
    input:
        counts=expand(join("output", "learn", "kmer-counts-{nb}.csv"), nb=FAS),
        base_counts=expand("{bf}", bf=base_counts),
    output:
        totals=join("output", "learn", "kmer-counts-total.csv"),
    log:
        join(out_dir, "learn", "log", "merge.log"),
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
        counts_files (list): List of paths to the CSV files containing k-mer counts.
        base_counts_path (str): Path to the base CSV file for merging.
        output_path (str): Path to save the merged dataframe.
        running_merge (DataFrame or None): Merged dataframe from counts_files.
        base_check (bool): Flag to check if the base file exists and can be merged.
        base_df (DataFrame or None): Dataframe loaded from base_counts_path.
    """

            def __init__(self, counts_files, base_counts_path, output_path):
                self.counts_files = counts_files
                self.base_counts_path = base_counts_path
                self.output_path = output_path
                self.running_merge = None
                self.base_check = False
                self.base_df = None

            def merge_dataframes(self):
                """
        Merges dataframes from the list of kmer count files.

        Loads each dataframe from counts_files, then successively merges them
        into a running merged dataframe.
        """
                for file_num, f in enumerate(self.counts_files):
                    kmer_counts = pd.read_csv(
                        f,
                        index_col="__index_level_0__",
                        header=0,
                        engine="pyarrow",
                        na_values=[""],
                    )
                    kmer_counts.fillna(0, inplace=True)
                    if file_num == 0:
                        self.running_merge = kmer_counts
                    else:
                        self.running_merge = (
                            pd.concat([self.running_merge, kmer_counts])
                            .reset_index()
                            .groupby("__index_level_0__", sort=False)
                            .sum(min_count=1)
                        ).fillna(0)
                    print(
                        f"Dataframes merged: {file_num} out of {len(self.counts_files)}"
                    )

            def check_for_base_file(self):
                """
        Checks for the presence of a base file to merge with.

        Sets the base_check flag to True if a CSV file is detected in the base path.
        """
                print("\nChecking for base file to merge with.\n")
                if "csv" in str(self.base_counts_path):
                    print(
                        "CSV detected. Matching annotations, kmers, and totals will be summed. New annotations and kmers will be added."
                    )
                    self.base_check = True
                elif self.base_counts_path == "":
                    print("No base directory detected\n")
                elif str(self.base_counts_path) == "input/base":
                    print("Empty base directory detected\n")
                else:
                    print(
                        "No file type detected. Please use a .csv file in input/base directory.\n"
                    )

            def confirm_kmer_counts_and_alphabet(self):
                """
        Confirms consistency between the alphabets and k-mer lengths
        of the running merged dataframe and the base dataframe.

        If any inconsistency is found, the base_check flag is set to False.
        """
                if self.base_check:
                    self.base_df = pd.read_csv(
                        str(self.base_counts_path),
                        index_col="__index_level_0__",
                        header=0,
                        engine="pyarrow",
                    )
                    print("\nBase Database: \n")
                    print(self.base_df)
                    check_1 = len(self.running_merge.columns.values)
                    alphabet_initial = set(
                        itertools.chain(
                            *[
                                list(x)
                                for x in self.running_merge.columns.values[3:check_1]
                            ]
                        )
                    )
                    alphabet_base = set(
                        itertools.chain(
                            *[list(x) for x in self.base_df.columns.values[3:check_1]]
                        )
                    )
                    if alphabet_base != alphabet_initial:
                        self.base_check = False
                        print("Different Alphabets Detected. Base File not merged.")
                    if len(str(self.running_merge.columns.values[1])) != len(
                        str(self.base_df.columns.values[1])
                    ):
                        self.base_check = False
                        print("Different kmer lengths detected. Base File not merged.")

            def merge_with_base(self):
                """
        Merges the running merged dataframe with the base dataframe,
        if the base_check flag is True.

        If the flag is False, only the running merged dataframe is saved to output.
        """
                if self.base_check:
                    print("\nMerged Database \n")
                    xy = (
                        pd.concat([self.base_df, self.running_merge])
                        .reset_index()
                        .groupby("__index_level_0__", sort=False)
                        .sum(min_count=1)
                    ).fillna(0)
                    xy_out = pa.Table.from_pandas(xy, preserve_index=True)
                    csv.write_csv(xy_out, self.output_path)
                    print(xy)
                else:
                    print("\nDatabase Merged. Not merged with base file.\n")
                    running_merge_out = pa.Table.from_pandas(
                        self.running_merge, preserve_index=True
                    )
                    csv.write_csv(running_merge_out, self.output_path)

            def execute_all(self):
                """
        Executes all the merging steps in sequence.

        This includes:
            1. Merging individual count dataframes.
            2. Checking for a base file.
            3. Confirming kmer counts and alphabet consistency.
            4. Merging with the base file if applicable.
        """
                self.merge_dataframes()
                self.check_for_base_file()
                self.confirm_kmer_counts_and_alphabet()
                self.merge_with_base()


        merger = Merge(input.counts, input.base_counts, output.totals)
        merger.execute_all()


rule eval_apply:
    input:
        data=join(
            "output",
            "vector"
            if config["learnapp"]["fragmentation"] == False
            else "vector_frag",
            "{nb}.npz",
        ),
        annotation=expand("{an}", an=annot_files),
        compare_associations=join("output", "learn", "kmer-counts-total.csv"),
    output:
        apply=join(
            "output",
            "eval_apply"
            if config["learnapp"]["fragmentation"] == False
            else "eval_apply_frag",
            "seq-annotation-scores-{nb}.csv",
        ),
    log:
        join(out_dir, "eval_apply", "log", "{nb}.log"),
    run:
        start_time = datetime.now()
        with open(log[0], "a") as f:
            f.write(f"start time:\t{start_time}\n")


        class KmerCompare:
            """
    Initializes the KmerCompare object.

    This object is designed to compare kmer counts with provided annotations.

    Attributes:
        compare_associations (str): Path to a CSV file containing kmer counts totals matrix.
        annotation_files (list): List of paths to files containing sequence annotations.
        input_data (str): Path to input data for kmer analysis.
        output_path (str): Path to save the result.
        annotation (list): List of dataframes loaded from annotation_files.
        kmer_count_totals (DataFrame or None): DataFrame of kmer count totals from compare_associations.
        seq_annot (dict): Dictionary mapping sequence IDs to annotations.
        kmerlist (list): List of unique kmers found in input data.
        seq_kmer_dict (dict): Dictionary mapping sequence IDs to their kmer counts.
        total_seqs (int): Total number of sequences processed.
        kmer_totals (list): Total counts for each kmer across all sequences.
    """

            def __init__(
                self, compare_associations, annotation_files, input_data, output_path
            ):
                self.compare_associations = compare_associations
                self.annotation_files = annotation_files
                self.input_data = input_data
                self.output_path = output_path
                self.annotation = []
                self.kmer_count_totals = None
                self.seq_annot = {}
                self.kmerlist = []
                self.seq_kmer_dict = {}
                self.total_seqs = 0
                self.kmer_totals = []

            def generate_inputs(self):
                """
        Generates the necessary inputs for comparison.

        Loads kmer counts and annotations into appropriate data structures.
        """
                self.kmer_count_totals = pd.read_csv(
                    str(self.compare_associations),
                    index_col="__index_level_0__",
                    header=0,
                    engine="c",
                )
                for f in self.annotation_files:
                    self.annotation.append(pd.read_table(f))
                seqs = self.annotation[0]["id"].tolist()
                anns = self.annotation[0]["TIGRFAMs"].tolist()
                for i, seqid in enumerate(seqs):
                    self.seq_annot[seqid] = anns[i]

            def generate_kmer_counts(self):
                """
        Generates a dictionary of kmer counts for each sequence.

        Processes the input data to count the occurrence of each kmer in each sequence.
        """
                kmerlist, df = skm.io.load_npz(self.input_data)
                self.kmerlist = kmerlist[0]
                seqids = df["sequence_id"]
                self.kmer_totals = [0] * len(self.kmerlist)
                k_len = len(self.kmerlist[0])

                for i, seq in enumerate(seqids):
                    v = df["sequence"][i]
                    k_counts = {}
                    items = [
                        v[item : (item + k_len)]
                        for item in range(0, (len((v)) - k_len + 1))
                    ]
                    for j in items:
                        k_counts[j] = k_counts.get(j, 0) + 1
                    store = [
                        k_counts[item] if item in k_counts else 0
                        for item in self.kmerlist
                    ]
                    for i, item in enumerate(self.kmerlist):
                        if item in k_counts:
                            self.kmer_totals[i] += k_counts[item]
                    self.seq_kmer_dict[seq] = store

            def add_known_unknown_tag(self):
                """
        Modifies sequence IDs to include a known or unknown tag based on annotation.

        Uses annotations to determine whether a sequence is known or unknown.
        Appends a tag to the sequence ID accordingly.
        """
                seqs = set(self.annotation[0]["id"].tolist())
                count = 0
                for seqid in list(self.seq_kmer_dict):
                    x = re.findall(r"\|(.*?)\|", seqid)[0]
                    if x not in seqs:
                        self.seq_kmer_dict[
                            (x + "_unknown_" + str(count))
                        ] = self.seq_kmer_dict.pop(seqid)
                    else:
                        self.seq_kmer_dict[
                            (self.seq_annot[x] + "_known_" + str(count))
                        ] = self.seq_kmer_dict.pop(seqid)
                    count += 1
                self.total_seqs = len(self.seq_kmer_dict)

            def construct_kmer_counts_dataframe(self):
                """
        Constructs a pandas DataFrame of kmer counts for each sequence.

        Returns:
            DataFrame: A DataFrame where rows represent sequences (and a total row),
                    and columns represent kmers.
        """
                kmer_counts = pd.DataFrame(self.seq_kmer_dict.values())
                kmer_counts.insert(0, "Annotations", 1, True)
                self.kmer_totals.insert(0, self.total_seqs)
                kmer_counts = pd.DataFrame(
                    np.insert(kmer_counts.values, 0, values=self.kmer_totals, axis=0)
                )
                kmer_counts.columns = ["Sequence count"] + list(self.kmerlist)
                kmer_counts.index = ["Totals"] + list(self.seq_kmer_dict.keys())
                return kmer_counts

            def match_kmer_counts_format(self, kmer_counts):
                """
        Matches the format of the provided kmer counts DataFrame to the format of the
        comparison data. Ensures columns align correctly.

        Args:
            kmer_counts (DataFrame): DataFrame of kmer counts to format.

        Returns:
            DataFrame: Formatted kmer counts DataFrame.
        """
                if len(str(kmer_counts.columns.values[10])) == len(
                    str(self.kmer_count_totals.columns.values[10])
                ):
                    compare_check = True
                else:
                    compare_check = False

                if compare_check:
                    check_1 = len(kmer_counts.columns.values)
                    alphabet_initial = set(
                        itertools.chain(
                            *[list(x) for x in kmer_counts.columns.values[10:check_1]]
                        )
                    )
                    alphabet_compare = set(
                        itertools.chain(
                            *[
                                list(x)
                                for x in self.kmer_count_totals.columns.values[
                                    10:check_1
                                ]
                            ]
                        )
                    )
                    if alphabet_compare != alphabet_initial:
                        compare_check = False

                if not compare_check:
                    print("Compare Check Failed. ")
                    sys.exit()

                kmer_counts.drop("Totals", axis=0, inplace=True)
                kmer_counts.drop("Sequence count", axis=1, inplace=True)

                self.kmer_count_totals.drop("Totals", axis=0, inplace=True)
                self.kmer_count_totals.drop("Kmer Count", axis=1, inplace=True)
                self.kmer_count_totals.drop("Sequence count", axis=1, inplace=True)

                column_order = list(
                    set(kmer_counts.columns) | set(self.kmer_count_totals.columns)
                )
                kmer_counts = kmer_counts.reindex(columns=column_order, fill_value=0)
                self.kmer_count_totals = self.kmer_count_totals.reindex(
                    columns=column_order, fill_value=0
                )

                return kmer_counts

            def calculate_cosine_similarity(self, kmer_counts):
                """
        Calculates the cosine similarity between the input kmer counts and comparison data.

        Args:
            kmer_counts (DataFrame): DataFrame of kmer counts for comparison.

        Returns:
            DataFrame: A DataFrame of cosine similarity scores.
        """
                cosine_df = sklearn.metrics.pairwise.cosine_similarity(
                    self.kmer_count_totals, kmer_counts
                ).T
                final_matrix_with_scores = pd.DataFrame(
                    cosine_df,
                    columns=self.kmer_count_totals.index,
                    index=kmer_counts.index,
                )
                return final_matrix_with_scores

            def filter_top_two_values(self, final_matrix_with_scores):
                """
        Filters the similarity scores to keep only the top two values for each row.

        Args:
            final_matrix_with_scores (DataFrame): DataFrame of similarity scores.

        Returns:
            DataFrame: DataFrame with all but the top two scores set to NaN.
        """
                top_2_indices = np.argsort(-final_matrix_with_scores.values, axis=1)[
                    :, :2
                ]
                mask = np.zeros_like(final_matrix_with_scores.values, dtype=bool)
                for i, (index_1, index_2) in enumerate(top_2_indices):
                    mask[i, index_1] = True
                    mask[i, index_2] = True
                final_matrix_with_scores.values[~mask] = np.nan
                return final_matrix_with_scores

            def write_output(self, final_matrix_with_scores):
                """
        Writes the provided DataFrame to a CSV file at the specified output path.

        Args:
            final_matrix_with_scores (DataFrame): DataFrame to write to CSV.
        """
                final_matrix_with_scores_write = pa.Table.from_pandas(
                    final_matrix_with_scores
                )
                csv.write_csv(final_matrix_with_scores_write, self.output_path)

            def execute_all(self, config):
                """
        Executes all the comparison steps in sequence.

        This includes:
            1. Loading inputs.
            2. Generating kmer counts.
            3. Adding known/unknown tags.
            4. Constructing a k-mer counts dataframe.
            5. Matching format with comparison data.
            6. Calculating cosine similarity.
            7. Filtering to keep top two values (if applicable).
            8. Writing results to output.
        """
                self.generate_inputs()
                self.generate_kmer_counts()
                self.add_known_unknown_tag()
                kmer_counts = self.construct_kmer_counts_dataframe()
                kmer_counts = self.match_kmer_counts_format(kmer_counts)
                final_matrix_with_scores = self.calculate_cosine_similarity(kmer_counts)
                if not config["learnapp"]["save_apply_associations"]:
                    final_matrix_with_scores = self.filter_top_two_values(
                        final_matrix_with_scores
                    )
                self.write_output(final_matrix_with_scores)


        analysis = KmerCompare(
            input.compare_associations, input.annotation, input.data, output.apply
        )
        analysis.execute_all(config)
        skm.utils.log_runtime(log[0], start_time)


rule evaluate:
    input:
        eval_apply_data=expand(
            join(
                "output",
                "eval_apply"
                if config["learnapp"]["fragmentation"] == False
                else "eval_apply_frag",
                "seq-annotation-scores-{nb}.csv",
            ),
            nb=FAS,
        ),
        base_confidence=expand("{bc}", bc=base_confidence),
    output:
        eval_conf="output/eval_conf/confidence-matrix.csv",
        eval_glob="output/eval_conf/global-confidence-scores.csv",
    params:
        modifer=config["learnapp"]["conf_weight_modifier"],
    log:
        join(out_dir, "eval_conf", "log", "conf.log"),
    run:
        start_time = datetime.now()
        with open(log[0], "a") as f:
            f.write(f"start time:\t{start_time}\n")


        class Evaluator:
            """
    The Evaluater class processes compares predictions with results and generates confidence metrics.

    Attributes:
        input_data (list): List of paths to input files.
        output_conf_path (str): Path to save the output confidence table.
        output_glob_path (str): Path to save the output global crosstab.
        confidence_data (list, optional): List of paths to base confidence files. Defaults to None.
    """

            def __init__(
                self,
                input_data,
                output_conf_path,
                output_glob_path,
                confidence_data=None,
            ):
                """
        Initializes the Evaluator object with input data and paths.
        """
                self.input_data = input_data
                self.confidence_data = confidence_data
                self.output_conf = output_conf_path
                self.output_glob = output_glob_path
                self.true_running_crosstab = None
                self.false_running_crosstab = None

            def read_and_transform_input_data(self, file_path):
                """
        Reads and transforms input data from a given CSV file.

        Args:
            file_path (str): Path to the input CSV file.

        Returns:
            tuple: Parsed and transformed values from the input data.
        """
                seq_ann_scores = pd.read_csv(
                    file_path,
                    index_col="__index_level_0__",
                    header=0,
                    engine="c",
                )
                max_value_index = seq_ann_scores.idxmax(axis="columns")
                result = max_value_index.keys()
                tf, known = self._generate_tf_and_known(max_value_index, result)

                seq_ann_vals = seq_ann_scores.values[
                    np.arange(len(seq_ann_scores))[:, None],
                    np.argpartition(-seq_ann_scores.values, np.arange(2), axis=1)[
                        :, :2
                    ],
                ]
                return seq_ann_vals, max_value_index, result, tf, known

            def _generate_tf_and_known(self, max_value_index, result):
                """
        Helper function to generate "True/False" and "Known/Unknown" labels based on input data.

        Args:
            max_value_index (Series): Maximum value's index from the data.
            result (Index): Result keys from data.

        Returns:
            tuple: Lists of "True/False" and "Known/Unknown" labels.
        """
                tf = []
                known = []
                for i, item in enumerate(list(max_value_index)):
                    if item in result[i]:
                        tf.append("T")
                    else:
                        tf.append("F")
                    if "unknown" in result[i]:
                        known.append("Unknown")
                    else:
                        known.append("Known")
                return tf, known

            def generate_diff_dataframe(
                self, seq_ann_vals, max_value_index, result, tf, known
            ):
                """
        Generates a dataframe showing the differences and other metrics.

        Args:
            seq_ann_vals (array-like): Array of sequence annotation values.
            max_value_index (Series): Maximum value's index from the data.
            result (Index): Result keys from data.
            tf (list): List of "True/False" labels.
            known (list): List of "Known/Unknown" labels.

        Returns:
            DataFrame: Constructed dataframe with computed differences and other metrics.
        """
                diff_df = pd.DataFrame(seq_ann_vals, columns=["Top", "Second"])
                diff_df["Difference"] = -(
                    np.diff(seq_ann_vals, axis=1).round(decimals=2)
                )
                diff_df["Prediction"] = list(max_value_index)
                diff_df["Actual"] = result
                diff_df["T/F"] = tf
                diff_df["Known/Unknown"] = known
                return diff_df

            def create_tf_crosstabs(self, diff_df):
                """
        Creates crosstabs for True and False predictions based on the given dataframe.

        Args:
            diff_df (DataFrame): DataFrame with differences and metrics.

        Returns:
            tuple: Crosstabs for True and False predictions.
        """
                known_true_diff_df = diff_df[
                    (diff_df["Known/Unknown"] == "Known") & (diff_df["T/F"] == "T")
                ]
                known_false_diff_df = diff_df[
                    (diff_df["Known/Unknown"] == "Known") & (diff_df["T/F"] == "F")
                ]
                possible_vals = [round(x * 0.01, 2) for x in range(0, 101)]
                true_crosstab = pd.crosstab(
                    known_true_diff_df.Prediction, known_true_diff_df.Difference
                )
                false_crosstab = pd.crosstab(
                    known_false_diff_df.Prediction, known_false_diff_df.Difference
                )
                return true_crosstab, false_crosstab

            def handle_running_crosstabs(
                self, true_crosstab, false_crosstab, iteration
            ):
                """
        Handles and updates running crosstabs over iterations.

        Args:
            true_crosstab (DataFrame): Crosstab for True predictions.
            false_crosstab (DataFrame): Crosstab for False predictions.
            iteration (int): Current iteration.
        """
                possible_vals = [round(x * 0.01, 2) for x in range(0, 101)]

                if iteration == 0:
                    self.true_running_crosstab = true_crosstab
                    self.false_running_crosstab = false_crosstab
                else:
                    self.true_running_crosstab = (
                        pd.concat([self.true_running_crosstab, true_crosstab])
                        .reset_index()
                        .groupby("Prediction", sort=False)
                        .sum(min_count=1)
                    ).fillna(0)
                    self.false_running_crosstab = (
                        pd.concat([self.false_running_crosstab, false_crosstab])
                        .reset_index()
                        .groupby("Prediction", sort=False)
                        .sum(min_count=1)
                    ).fillna(0)

                add_to_true_df = pd.DataFrame(
                    0,
                    index=sorted(
                        set(self.false_running_crosstab.index)
                        - set(self.true_running_crosstab.index)
                    ),
                    columns=self.true_running_crosstab.columns,
                )
                add_to_false_df = pd.DataFrame(
                    0,
                    index=sorted(
                        set(self.true_running_crosstab.index)
                        - set(self.false_running_crosstab.index)
                    ),
                    columns=self.false_running_crosstab.columns,
                )

                self.true_running_crosstab = pd.concat(
                    [self.true_running_crosstab, add_to_true_df]
                )[
                    sorted(
                        list(
                            set(possible_vals) & set(self.true_running_crosstab.columns)
                        )
                    )
                ].assign(
                    **dict.fromkeys(
                        list(
                            map(
                                str,
                                sorted(
                                    list(
                                        set(possible_vals)
                                        ^ set(
                                            self.true_running_crosstab.columns.astype(
                                                float
                                            )
                                        )
                                    )
                                ),
                            )
                        ),
                        0,
                    )
                )
                self.false_running_crosstab = pd.concat(
                    [self.false_running_crosstab, add_to_false_df]
                )[
                    sorted(
                        list(
                            set(possible_vals)
                            & set(self.false_running_crosstab.columns)
                        )
                    )
                ].assign(
                    **dict.fromkeys(
                        list(
                            map(
                                str,
                                sorted(
                                    list(
                                        set(possible_vals)
                                        ^ set(
                                            self.false_running_crosstab.columns.astype(
                                                float
                                            )
                                        )
                                    )
                                ),
                            )
                        ),
                        0,
                    )
                )

                self.true_running_crosstab.index.names = ["Prediction"]
                self.false_running_crosstab.index.names = ["Prediction"]
                self.true_running_crosstab.sort_index(inplace=True)
                self.false_running_crosstab.sort_index(inplace=True)
                self.true_running_crosstab.columns = (
                    self.true_running_crosstab.columns.astype(float)
                )
                self.false_running_crosstab.columns = (
                    self.false_running_crosstab.columns.astype(float)
                )
                self.true_running_crosstab = self.true_running_crosstab[
                    sorted(self.true_running_crosstab.columns)
                ]
                self.false_running_crosstab = self.false_running_crosstab[
                    sorted(self.false_running_crosstab.columns)
                ]

            def generate_inputs(self):
                """
        Generates and processes input data for each file in the input_data attribute.
        """
                for j, f in enumerate(self.input_data):
                    (
                        seq_ann_vals,
                        max_value_index,
                        result,
                        tf,
                        known,
                    ) = self.read_and_transform_input_data(f)
                    diff_df = self.generate_diff_dataframe(
                        seq_ann_vals, max_value_index, result, tf, known
                    )
                    true_crosstab, false_crosstab = self.create_tf_crosstabs(diff_df)
                    known_count = (
                        diff_df["Known/Unknown"].value_counts().get("Known", 0)
                    )
                    self.handle_running_crosstabs(true_crosstab, false_crosstab, j)
                return

            def generate_global_crosstab(self):
                """
        Generates a global crosstab based on the running True and False crosstabs.

        Returns:
            DataFrame: Computed global crosstab.
        """
                return self.true_running_crosstab / (
                    self.true_running_crosstab + self.false_running_crosstab
                )

            def calculate_distributions(self):
                """
        Calculates distributions for True and False predictions.

        Returns:
            tuple: Total distributions for True and False predictions.
        """
                true_total_dist = self.true_running_crosstab.sum(
                    numeric_only=True, axis=0
                )
                false_total_dist = self.false_running_crosstab.sum(
                    numeric_only=True, axis=0
                )

                return true_total_dist, false_total_dist

            def compute_ratio_distribution(self, true_total_dist, false_total_dist):
                """
        Computes the ratio distribution for True and False total distributions.

        Args:
            true_total_dist (Series): Total distribution for True predictions.
            false_total_dist (Series): Total distribution for False predictions.

        Returns:
            Series: Computed ratio distribution.
        """
                ratio_total_dist = true_total_dist / (
                    true_total_dist + false_total_dist
                )
                return ratio_total_dist.interpolate(method="linear")

            def check_confidence_merge(self, new_ratio_dist):
                """
        Merge with base confidence file if available.

        Args:
            new_ratio_dist (Series): Total distribution for T/(T+F) predictions.

        Returns:
            DataFrame: Updated ratio distribution.
        """
                sum_series = (
                    self.true_running_crosstab.sum() + self.false_running_crosstab.sum()
                )
                current_weight = (
                    self.true_running_crosstab.values.sum()
                    + self.false_running_crosstab.values.sum()
                )

                updated_data = pd.DataFrame(new_ratio_dist, columns=["confidence"])
                updated_data = updated_data.assign(
                    weight=[current_weight] * 101, sum=sum_series
                )

                if self.confidence_data and len(self.confidence_data) == 1:
                    prior_conf = pd.read_csv(
                        self.confidence_data[0], index_col="Difference"
                    )
                    print(f"Prior Confidence Data:\n{prior_conf}")

                    total_weight_prior = prior_conf["weight"]
                    k_factor = 1 + params.modifer * (
                        current_weight / (current_weight + total_weight_prior)
                    )
                    out_weight = total_weight_prior + current_weight
                    weighted_current = k_factor * current_weight
                    total_weight = total_weight_prior + weighted_current
                    prior_weighted_score = (
                        prior_conf["confidence"] * prior_conf["weight"]
                    )
                    current_weighted_score = (
                        updated_data["confidence"] * weighted_current
                    )

                    updated_confidence = (
                        prior_weighted_score + current_weighted_score
                    ) / total_weight

                    updated_data = pd.DataFrame({"confidence": updated_confidence})

                    updated_data = updated_data.assign(
                        weight=out_weight,
                        sum=sum_series + prior_conf["sum"],
                    )

                    print(f"Final Confidence Data\n{updated_data}")
                else:
                    print(
                        "Base confidence file not found. Only one file is allowed in base/confidence."
                    )

                return updated_data

            def save_results(
                self,
                ratio_total_dist,
                true_total_dist,
                false_total_dist,
                output_path,
                conf_path,
            ):
                """
        Saves the computed results to specified paths.

        Args:
            ratio_total_dist (Series): Ratio distribution to be saved.
            output_path (str): Path to save the ratio distribution.
            conf_path (str): Path to save the confidence table.
        """

                ratio_total_dist.to_csv(output_path)
                table = pa.Table.from_pandas(self.generate_global_crosstab())
                csv.write_csv(table, conf_path)

            def execute_all(self):
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
                ratio_crosstab = self.generate_global_crosstab()
                true_dist, false_dist = self.calculate_distributions()
                ratio_dist = self.compute_ratio_distribution(true_dist, false_dist)
                ratio_dist = self.check_confidence_merge(ratio_dist)
                self.save_results(
                    ratio_dist,
                    true_dist,
                    false_dist,
                    self.output_glob,
                    self.output_conf,
                )


        evaluator = Evaluator(
            input.eval_apply_data,
            output.eval_conf,
            output.eval_glob,
            input.base_confidence,
        )
        evaluator.execute_all()


        skm.utils.log_runtime(log[0], start_time)
