# ---------------------------------------------------------
# Imports
# ---------------------------------------------------------

import re
from os.path import basename, splitext, join

import numpy as np
import pandas as pd

from typing import Dict, List, Optional

import snekmer as skm
from snekmer.learn import generate_kmer_counts


# ---------------------------------------------------------
# Run logic
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
    out_dir (str): Base output directory where 'learn/' subdir will be used.
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
    out_dir: Optional[str]

    # --- Optimization: compile regex once at class level ---
    _SEQID_PATTERN = re.compile(r"\|(.*?)\|")

    def __init__(self, out_dir: Optional[str] = None):
        self.annotation = []
        self.seq_annot = {}
        self.kmer_list = []
        self.df = None
        self.seq_ids = []
        self.kmer_totals = []
        self.seq_kmer_dict = {}
        self.annotation_counts = {}
        self.total_seqs = 0
        self.out_dir = out_dir

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

        self.seq_annot = dict(zip(seqs, anns))
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
        self.kmer_totals = [0] * len(self.kmer_list)

    def filter_and_construct(self) -> None:
        """
        Filters sequences not present in annotations and constructs annotation counts.
        """
        self.total_seqs = len(self.seq_kmer_dict)
        pattern = self._SEQID_PATTERN

        # --- Optimization: convert values to numpy arrays for fast aggregation,
        #     build result dict directly instead of repeated pop/reassign ---
        aggregated: Dict[str, np.ndarray] = {}
        self.annotation_counts = {}

        for seqid, kmer_vals in self.seq_kmer_dict.items():
            x = pattern.findall(seqid)[0]
            if x not in self.seqs:
                continue

            ann = self.seq_annot[x]
            arr = np.asarray(kmer_vals)

            if ann not in aggregated:
                aggregated[ann] = arr
                self.annotation_counts[ann] = 1
            else:
                aggregated[ann] += arr
                self.annotation_counts[ann] += 1

        # Convert back to lists to preserve downstream DataFrame construction behavior
        self.seq_kmer_dict = {k: v.tolist() for k, v in aggregated.items()}

    def format_and_write_output(self, input_data: str) -> None:
        """
        Writes processed kmer counts to an output CSV file.

        Args:
            input_data (str): Path to the data file (used for naming the output file).
        """
        if self.out_dir is None:
            raise ValueError(
                "Library.out_dir is not set. Please initialize Library(out_dir=...) "
                "or use run_learn(...) which sets it for you."
            )

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
        base = basename(input_data)
        name = splitext(base)[0]
        out_name = join(self.out_dir, "learn", f"kmer_counts_{name}.csv")
        kmer_counts.index.name = "__index_level_0__"
        kmer_counts.to_csv(out_name, index=True)

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
        self.kmer_list, self.kmer_totals, self.seq_kmer_dict = generate_kmer_counts(
            input_data,
            self.kmer_list,
            self.kmer_totals,
            self.seq_kmer_dict,
            False,
        )
        self.filter_and_construct()
        self.format_and_write_output(input_data)


# ---------------------------------------------------------
# Public helper for demos / external calls
# ---------------------------------------------------------

def run_learn(
    input_annotation: List[str],
    input_data: str,
    config: Dict,
    out_dir: Optional[str] = None,
) -> None:
    """
    Run the Learn step outside of Snakemake.

    Parameters
    ----------
    input_annotation : list[str]
        Annotation files (TSV) with at least columns 'id' and 'family'.
    input_data : str
        Path to the NPZ file containing k-merized sequences.
    config : dict
        Configuration dict with keys 'alphabet', 'k', and 'nested_output'
        (same as Snakemake config["alphabet"], config["k"], config["nested_output"]).
    out_dir : str, optional
        Base output directory. If None, it will be computed via
        `snekmer.io.define_output_dir(config["alphabet"], config["k"], nested=config["nested_output"])`.
    """
    if out_dir is None:
        out_dir = skm.io.define_output_dir(
            config["alphabet"], config["k"], nested=config["nested_output"]
        )

    library = Library(out_dir=out_dir)
    library.execute_all(input_annotation, input_data)


# ---------------------------------------------------------
# Snakemake entry point
# ---------------------------------------------------------

# This block keeps the existing Snakemake behavior intact.
# Snakemake injects a global `snakemake` object when using `script:`.
if "snakemake" in globals():
    config = snakemake.config
    out_dir = skm.io.define_output_dir(
        config["alphabet"], config["k"], nested=config["nested_output"]
    )

    library = Library(out_dir=out_dir)
    library.execute_all(snakemake.input.annotation, snakemake.input.data)