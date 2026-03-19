# ---------------------------------------------------------
# Imports
# ---------------------------------------------------------

import itertools

import pandas as pd
import pyarrow as pa
import pyarrow.csv as csv

import snekmer as skm


# ---------------------------------------------------------
# Core logic
# ---------------------------------------------------------

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
    base_kmer_counts (DataFrame or None): Dataframe loaded from base_counts_path.
    """
    counts_files: list
    base_counts_path: str
    output_path: str
    running_merge: pd.DataFrame
    base_check: bool
    base_kmer_counts: None

    def __init__(self, counts_files: list, base_counts_path: str, output_path: str) -> None:
        self.counts_files = counts_files
        self.base_counts_path = base_counts_path
        self.output_path = output_path
        self.running_merge = None
        self.base_check = False
        self.base_kmer_counts = None

    def merge_dataframes(self) -> None:
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
                f"Dataframes merged: {file_num + 1} out of {len(self.counts_files)}"
            )

    def check_for_base_file(self) -> None:
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

    def confirm_kmer_counts_and_alphabet(self) -> None:
        """
        Confirms consistency between the alphabets and k-mer lengths
        of the running merged dataframe and the base dataframe.

        If any inconsistency is found, the base_check flag is set to False.
        """
        if self.base_check:
            self.base_kmer_counts = pd.read_csv(
                str(self.base_counts_path),
                index_col="__index_level_0__",
                header=0,
                engine="pyarrow",
            )
            print("\nBase Database: \n")
            print(self.base_kmer_counts)
            check = len(self.running_merge.columns.values)
            alphabet_initial = set(
                itertools.chain(
                    *[
                        list(x)
                        for x in self.running_merge.columns.values[3:check]
                    ]
                )
            )
            alphabet_base = set(
                itertools.chain(
                    *[
                        list(x)
                        for x in self.base_kmer_counts.columns.values[3:check]
                    ]
                )
            )
            if alphabet_base != alphabet_initial:
                self.base_check = False
                print("Different Alphabets Detected. Base File not merged.")
            if len(str(self.running_merge.columns.values[1])) != len(
                str(self.base_kmer_counts.columns.values[1])
            ):
                self.base_check = False
                print("Different kmer lengths detected. Base File not merged.")

    def merge_with_base(self) -> None:
        """
        Merges the running merged dataframe with the base dataframe,
        if the base_check flag is True.

        If the flag is False, only the running merged dataframe is saved to output.
        """
        if self.base_check:
            print("\nMerged Database \n")
            data = (
                pd.concat([self.base_kmer_counts, self.running_merge])
                .reset_index()
                .groupby("__index_level_0__", sort=False)
                .sum(min_count=1)
            ).fillna(0)
            file_out = pa.Table.from_pandas(data, preserve_index=True)
            csv.write_csv(file_out, self.output_path)
        else:
            print("\nDatabase Merged. Not merged with base file.\n")
            running_merge_out = pa.Table.from_pandas(
                self.running_merge, preserve_index=True
            )
            csv.write_csv(running_merge_out, self.output_path)

    def execute_all(self) -> None:
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


# ---------------------------------------------------------
# Public helper for demos / external calls
# ---------------------------------------------------------

def run_merge(counts_files, base_counts_path, output_path) -> None:
    """
    Run the Merge step outside of Snakemake.

    Parameters
    ----------
    counts_files : list[str]
        Paths to k-mer count CSV files to be merged.
    base_counts_path : str
        Path to an existing base CSV (or empty string / dummy path if none).
    output_path : str
        Output CSV path for the merged totals.
    """
    merger = Merge(counts_files, base_counts_path, output_path)
    merger.execute_all()


# ---------------------------------------------------------
# Snakemake entry point
# ---------------------------------------------------------

# Keep Snakemake behavior unchanged; only run this block when Snakemake injects
# the global `snakemake` object.
if "snakemake" in globals():
    merger = Merge(
        snakemake.input.counts,
        snakemake.input.base_counts,
        snakemake.output.totals,
    )
    merger.execute_all()
