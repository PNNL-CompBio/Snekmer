# ---------------------------------------------------------
# Imports
# ---------------------------------------------------------

import gzip
from typing import List, Dict, Any

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.csv as csv
from sklearn.metrics.pairwise import cosine_similarity

import snekmer as skm
from snekmer.learn import generate_kmer_counts, match_kmer_counts_format


# ---------------------------------------------------------
# Core class
# ---------------------------------------------------------


class CompareReverseSeqs:
    """
    Compare reverse-sequence kmer counts against a learned kmer-association matrix.

    Attributes
    ----------
    compare_associations : str
        Path to a CSV file containing the kmer totals/associations matrix.
    annotation_files : list of str
        Paths to annotation files (TSV) – currently only loaded to mirror the
        forward-sequence pipeline; not strictly required here for scoring.
    input_data : str
        Path to the NPZ (or equivalent) with sequence + kmerized data.
    output_path : str
        Path to write the gzipped CSV of cosine similarity scores.
    """

    compare_associations: str
    annotation_files: List[str]
    input_data: str
    output_path: str
    annotation: List[pd.DataFrame]
    kmer_count_totals: pd.DataFrame
    seq_annot: Dict[str, Any]
    kmer_list: List[str]
    seq_kmer_dict: Dict[str, list]
    total_seqs: int
    kmer_totals: List[int]

    def __init__(
        self,
        compare_associations: str,
        annotation_files: List[str],
        input_data: str,
        output_path: str,
    ) -> None:
        self.compare_associations = compare_associations
        self.annotation_files = annotation_files
        self.input_data = input_data
        self.output_path = output_path

        self.annotation = []
        self.kmer_count_totals = None
        self.seq_annot = {}
        self.kmer_list = []
        self.seq_kmer_dict = {}
        self.total_seqs = 0
        self.kmer_totals = []

    def generate_inputs(self) -> None:
        """
        Load association matrix and (optionally) annotation files.
        """
        self.kmer_count_totals = pd.read_csv(
            str(self.compare_associations),
            index_col="__index_level_0__",
            header=0,
            engine="c",
        )

        # Annotations currently not used downstream for reverse seqs, but we keep
        # behavior compatible with forward pipeline and potential future checks.
        for f in self.annotation_files:
            self.annotation.append(pd.read_table(f))

    def construct_kmer_counts_dataframe(self) -> pd.DataFrame:
        """
        Constructs a pandas DataFrame of kmer counts for each sequence.

        Returns
        -------
        DataFrame
            Rows = sequences (+ Totals row), columns = kmers.
        """
        kmer_counts = pd.DataFrame(self.seq_kmer_dict.values())
        kmer_counts.insert(0, "Annotations", 1, True)
        self.kmer_totals.insert(0, self.total_seqs)
        kmer_counts = pd.DataFrame(
            np.insert(kmer_counts.values, 0, values=self.kmer_totals, axis=0)
        )
        kmer_counts.columns = ["Sequence count"] + list(self.kmer_list)
        kmer_counts.index = ["Totals"] + list(self.seq_kmer_dict.keys())
        return kmer_counts

    def calculate_cosine_similarity(self, kmer_counts: pd.DataFrame) -> pd.DataFrame:
        """
        Calculates cosine similarity between input kmer counts and the association matrix.

        Parameters
        ----------
        kmer_counts : DataFrame
            Sequence-by-kmer count matrix (with Totals row).

        Returns
        -------
        DataFrame
            Cosine similarity scores (rounded to 3 decimals).
        """
        cosine_dataframe = cosine_similarity(self.kmer_count_totals, kmer_counts).T
        final_matrix_with_scores = pd.DataFrame(
            cosine_dataframe,
            columns=self.kmer_count_totals.index,
            index=kmer_counts.index,
        )
        final_matrix_with_scores = final_matrix_with_scores.round(3)
        return final_matrix_with_scores

    def write_output(self, final_matrix_with_scores: pd.DataFrame) -> None:
        """
        Write similarity matrix to gzipped CSV at self.output_path.
        """
        final_matrix_with_scores_write = pa.Table.from_pandas(
            final_matrix_with_scores
        )
        with gzip.open(self.output_path, "wb") as gzipped_file:
            csv.write_csv(final_matrix_with_scores_write, gzipped_file)

    def execute_all(self) -> None:
        """
        Full pipeline for reverse-sequence evaluation:

          1. Load association matrix (+ optional annotations).
          2. Generate reverse-sequence kmer counts using generate_kmer_counts(..., reverse=True).
          3. Construct kmer count DataFrame.
          4. Align to association matrix with match_kmer_counts_format.
          5. Compute cosine similarity.
          6. Write gzipped CSV to self.output_path.
        """
        # 1: load inputs
        self.generate_inputs()

        # 2: generate kmer counts in reverse mode (True)
        self.kmer_list, self.kmer_totals, self.seq_kmer_dict = generate_kmer_counts(
            self.input_data,
            self.kmer_list,
            self.kmer_totals,
            self.seq_kmer_dict,
            True,  # reverse sequences
        )
        self.total_seqs = len(self.seq_kmer_dict)

        # 3: build kmer count matrix
        kmer_counts = self.construct_kmer_counts_dataframe()

        # 4: align with association db
        kmer_counts, self.kmer_count_totals = match_kmer_counts_format(
            kmer_counts, self.kmer_count_totals
        )

        # 5: cosine similarity
        final_matrix_with_scores = self.calculate_cosine_similarity(kmer_counts)

        # 6: write to disk
        self.write_output(final_matrix_with_scores)


# ---------------------------------------------------------
# Snakemake entry point
# ---------------------------------------------------------

if "snakemake" in globals():
    analysis = CompareReverseSeqs(
        compare_associations=snakemake.input.compare_associations,
        annotation_files=list(snakemake.input.annotation),
        input_data=snakemake.input.data,
        output_path=snakemake.output.apply,
    )
    analysis.execute_all()
