# ---------------------------------------------------------
# Imports
# ---------------------------------------------------------

import re
from os.path import join

import numpy as np
import pandas as pd

from typing import Dict, List, Optional

import snekmer as skm
from snekmer.learn import generate_kmer_counts


# ---------------------------------------------------------
# Files and Parameters
# ---------------------------------------------------------

config = snakemake.config


out_dir = skm.io.define_output_dir(
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
    seq_annot (dict): A dictionary mapping sequence IDs to their annotations.
    kmer_list (list): A list of unique kmers present in the data.
    df (DataFrame or None): DataFrame containing kmer data.
    seq_ids (list): A list of sequence IDs.
    kmer_totals (list): A list to store total counts of each k-mer across all sequences.
    seq_kmer_dict (dict): A dictionary mapping sequence IDs to their k-mer counts.
    annotation_counts (dict): A dictionary mapping annotations to their counts.
    total_seqs (int): The total number of sequences after filtering.
    """
    annotation: List[pd.DataFrame]
    seq_annot: Dict[str, str]
    kmer_list: List[str]
    df: Optional[pd.DataFrame]
    seq_ids: Optional[pd.Series]
    kmer_totals: List[int]
    seq_kmer_dict: Dict[str, int]
    annotation_counts: Dict[str, int]
    total_seqs: int

    def __init__(self):
        self.annotation = []
        self.seq_annot = {}
        self.kmer_list = []
        self.df = None
        self.seq_ids = []
        self.kmer_totals = []
        self.seq_kmer_dict = {}
        self.annotation_counts = {}
        self.total_seqs = 0

    def load_annotations(self, input_annotation: List[str]) -> None:
        """
        Load annotations from a list of provided input files.

        Args:
            input_annotation (list): List of file paths containing annotations.
        """
        for f in input_annotation:
            self.annotation.append(pd.read_table(f))
        annotations = pd.concat(self.annotation)
        seqs = annotations["id"].tolist()
        anns = [str(fam) for fam in annotations["family"].tolist()]

        for i, seqid in enumerate(seqs):
            self.seq_annot[seqid] = anns[i]
        self.seqs = set(seqs)

    def load_data(self, input_data: str) -> None:
        """
        Load and format kmer data from the provided input.

        Args:
            input_data (str): Path to the data file.
        """
        self.kmer_list, self.df = skm.io.load_npz(input_data)
        self.kmer_list = self.kmer_list[0]
        self.seq_ids = self.df["sequence_id"]
        for item in self.kmer_list:
            self.kmer_totals.append(0)

    def filter_and_construct(self) -> None:
        """
        Filters sequences not present in annotations and constructs annotation counts.
        """
        self.total_seqs = len(self.seq_kmer_dict)
        for i, seqid in enumerate(list(self.seq_kmer_dict)):
            x = re.findall(r"\|(.*?)\|", seqid)[
                0
            ]  # A Note, it could be useful to allow a user to define their own regex.
            if x not in self.seqs:
                del self.seq_kmer_dict[seqid]
            else:
                self.process_annotation_counts(seqid, x)

    def format_and_write_output(self, input_data: str) -> None:
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
        colnames = ["Sequence count"] + ["Kmer Count"] + list(self.kmer_list)
        kmer_counts = pd.DataFrame(
            np.insert(kmer_counts.values, 0, values=self.kmer_totals, axis=0)
        )
        kmer_counts.columns = colnames
        new_index = ["Totals"] + list(self.annotation_counts.keys())
        kmer_counts.index = new_index
        kmer_counts.replace(0, "", inplace=True)
        out_name = join(
            out_dir, "learn", "kmer-counts-" + str(input_data)[14:-4] + ".csv"
        )
        kmer_counts.index.name = "__index_level_0__"
        kmer_counts.to_csv(out_name, index=True)

    def process_annotation_counts(self, seqid: str, x: str) -> None:
        """
        Processes annotation counts by aggregating them based on annotation labels.

        Args:
            seqid (str): Sequence ID.
            x (str): Extracted annotation ID from seqid.
        """
        if self.seq_annot[x] not in self.seq_kmer_dict:
            self.seq_kmer_dict[self.seq_annot[x]] = self.seq_kmer_dict.pop(seqid)
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

    def execute_all(self, input_annotation: List[str], input_data: str) -> None:
        """
        Orchestrates the full k-mer counting and annotation workflow.

        This method runs the following steps in sequence:
          1. Loads annotation tables from the provided files.
          2. Loads k-mer data and associated sequence IDs.
          3. Computes k-mer counts for each sequence.
          4. Filters out unannotated sequences and aggregates counts by annotation.
          5. Formats the final k-mer count matrix and writes it to CSV.

        Args:
            input_annotation (List[str]): Paths to one or more annotation files (TSV).
            input_data (str): Path to the k-mer data file (e.g., NPZ).
        """
        self.load_annotations(input_annotation)
        self.load_data(input_data)
        self.kmer_list, self.kmer_totals, self.seq_kmer_dict = generate_kmer_counts(input_data, self.kmer_list, self.kmer_totals, self.seq_kmer_dict, False)
        self.filter_and_construct()
        self.format_and_write_output(input_data)


library = Library()
library.execute_all(snakemake.input.annotation, snakemake.input.data)
