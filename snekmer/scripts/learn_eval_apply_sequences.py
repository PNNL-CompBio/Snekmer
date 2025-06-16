# ---------------------------------------------------------
# Imports
# ---------------------------------------------------------

import gzip
import itertools
import re
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


class EvaluateSequences:
    """
    Initializes the EvaluateSequences object.

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
    kmer_count_totals: pd.DataFrame
    seq_annot: dict
    kmer_list: list
    seq_kmer_dict: dict
    total_seqs: int
    kmer_totals: list
    
    def __init__(
        self, compare_associations: str, annotation_files: list, input_data: str, output_path: str) -> None:
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
        anns = self.annotation[0]["Family"].tolist()

        for i, seqid in enumerate(seqs):
            self.seq_annot[seqid] = anns[i]

    def generate_kmer_counts(self) -> None:
        """
        Generates a dictionary of kmer counts for each sequence.

        Processes the input data to count the occurrence of each kmer in each sequence.
        """
        kmer_list, df = skm.io.load_npz(self.input_data)
        self.kmer_list = kmer_list[0]
        seq_ids = df["sequence_id"]
        self.kmer_totals = [0] * len(self.kmer_list)
        kmer_len = len(self.kmer_list[0])

        for i, seq in enumerate(seq_ids):
            v = df["sequence"][i]
            k_counts = {}
            items = [
                v[item : (item + kmer_len)]
                for item in range(0, (len((v)) - kmer_len + 1))
            ]
            for j in items:
                k_counts[j] = k_counts.get(j, 0) + 1
            store = [
                k_counts[item] if item in k_counts else 0
                for item in self.kmer_list
            ]
            for idx, item in enumerate(self.kmer_list):
                if item in k_counts:
                    self.kmer_totals[idx] += k_counts[item]
            self.seq_kmer_dict[seq] = store

    def add_known_tag(self) -> None:
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
                self.seq_kmer_dict[(str(x) + "_unknown_" + str(count))] = (
                    self.seq_kmer_dict.pop(seqid)
                )
            else:
                self.seq_kmer_dict[
                    (str(self.seq_annot[x]) + "_known_" + str(count))
                ] = self.seq_kmer_dict.pop(seqid)
            count += 1
        self.total_seqs = len(self.seq_kmer_dict)

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

    def calculate_cosine_similarity(self, kmer_counts: pd.DataFrame) -> None:
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
        return final_matrix_with_scores

    def filter_top_two_values(self, final_matrix_with_scores: pd.DataFrame) -> None:
        """
        Filters the similarity scores to keep only the top two values for each row.

        Args:
            final_matrix_with_scores (DataFrame): DataFrame of similarity scores.

        Returns:
            DataFrame: DataFrame with all but the top two scores set to NaN.
        """
        top_two_indices = np.argsort(-final_matrix_with_scores.values, axis=1)[:, :2]
        mask = np.zeros_like(final_matrix_with_scores.values, dtype=bool)
        for i, (index_one, index_two) in enumerate(top_two_indices):
            mask[i, index_one] = True
            mask[i, index_two] = True
        final_matrix_with_scores.values[~mask] = np.nan
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

    def execute_all(self, config: str) -> None:
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
        self.kmer_list, self.kmer_totals, self.seq_kmer_dict = generate_kmer_counts(self.input_data, 
                                                                                    self.kmer_list, 
                                                                                    self.kmer_totals, 
                                                                                    self.seq_kmer_dict, 
                                                                                    False)
        
        self.add_known_tag()
        kmer_counts = self.construct_kmer_counts_dataframe()
        kmer_counts, self.kmer_count_totals = match_kmer_counts_format(kmer_counts,self.kmer_count_totals)
        final_matrix_with_scores = self.calculate_cosine_similarity(kmer_counts)
        if not config["learnapp"]["save_apply_associations"]:
            final_matrix_with_scores = self.filter_top_two_values(
                final_matrix_with_scores
            )
        self.write_output(final_matrix_with_scores)


analysis = EvaluateSequences(
    snakemake.input.compare_associations,
    snakemake.input.annotation,
    snakemake.input.data,
    snakemake.output.apply
)

analysis.execute_all(config)
