# ---------------------------------------------------------
# Imports
# ---------------------------------------------------------

import gzip
import itertools
import sys

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.csv as csv
from sklearn.metrics.pairwise import cosine_similarity

import snekmer as skm
from snekmer.learn import generate_kmer_counts, match_kmer_counts_format

# ---------------------------------------------------------
# Files and Parameters
# ---------------------------------------------------------

config = snakemake.config

# ---------------------------------------------------------
# Run script
# ---------------------------------------------------------


class CompareReverseSeqs:
    """
    Initializes the CompareReverseSeqs object.

    This object is designed to compare kmer counts with provided annotations.

    Attributes:
    compare_associations (str): Path to a CSV file containing kmer counts totals matrix.
    annotation_files (list): List of paths to files containing sequence annotations.
    input_data (str): Path to input data for kmer analysis.
    output_path (str): Path to save the result.
    annotation (list): List of dataframes loaded from annotation_files.
    kmer_count_totals (DataFrame or None): DataFrame of kmer count totals from compare_associations.
    seq_annot (dict): Dictionary mapping sequence IDs to annotations.
    kmer_list (list): List of unique kmers found in input data.
    seq_kmer_dict (dict): Dictionary mapping sequence IDs to their kmer counts.
    total_seqs (int): Total number of sequences processed.
    kmer_totals (list): Total counts for each kmer across all sequences.
    """
    compare_associations: str
    annotation_files: list
    input_data: str
    output_path: str
    annotation: list
    kmer_count_totals: None
    seq_annot: dict
    kmer_list: list
    seq_kmer_dict: dict
    total_seqs: int
    kmer_totals: list
        
    def __init__(
        self, compare_associations: str, annotation_files: list, input_data: str, output_path: str
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
        anns = self.annotation[0]["family"].tolist()
        for i, seqid in enumerate(seqs):
            self.seq_annot[seqid] = anns[i]

    def construct_kmer_counts_dataframe(self) -> None:
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
        kmer_counts.columns = ["Sequence count"] + list(self.kmer_list)
        kmer_counts.index = ["Totals"] + list(self.seq_kmer_dict.keys())
        return kmer_counts

    def calculate_cosine_similarity(self, kmer_counts: pd.DataFrame) -> pd.DataFrame:
        """
        Calculates the cosine similarity between the input kmer counts and comparison data.

        Args:
            kmer_counts (DataFrame): DataFrame of kmer counts for comparison.

        Returns:
            DataFrame: A DataFrame of cosine similarity scores.
        """
        cosine_dataframe = cosine_similarity(
            self.kmer_count_totals, kmer_counts
        ).T
        final_matrix_with_scores = pd.DataFrame(
            cosine_dataframe,
            columns=self.kmer_count_totals.index,
            index=kmer_counts.index,
        )
        final_matrix_with_scores = final_matrix_with_scores.round(
            3
        )  # Adding rounding for storage space.
        return final_matrix_with_scores

    def write_output(self, final_matrix_with_scores: pd.DataFrame) -> None:
        """
        Writes the provided DataFrame to a CSV file at the specified output path.

        Args:
            final_matrix_with_scores (DataFrame): DataFrame to write to CSV.
        """
        final_matrix_with_scores_write = pa.Table.from_pandas(final_matrix_with_scores)
        with gzip.open(self.output_path, "wb") as gzipped_file:
            csv.write_csv(final_matrix_with_scores_write, gzipped_file)

    def execute_all(self) -> None:
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
        self.generate_inputs()
        self.kmer_list, self.kmer_totals, self.seq_kmer_dict = generate_kmer_counts(
            self.input_data, self.kmer_list, self.kmer_totals, self.seq_kmer_dict, True
        )

        kmer_counts = self.construct_kmer_counts_dataframe()
        kmer_counts, self.kmer_count_totals = match_kmer_counts_format(
            kmer_counts, self.kmer_count_totals
        )
        final_matrix_with_scores = self.calculate_cosine_similarity(kmer_counts)
        self.write_output(final_matrix_with_scores)


analysis = CompareReverseSeqs(
    snakemake.input.compare_associations,
    snakemake.input.annotation,
    snakemake.input.data,
    snakemake.output.apply
)
analysis.execute_all()
