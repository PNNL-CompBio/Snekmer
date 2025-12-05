# ---------------------------------------------------------
# Imports
# ---------------------------------------------------------

import gzip
import re
from typing import List, Optional, Dict, Any, Mapping

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.csv as csv
from sklearn.metrics.pairwise import cosine_similarity

import snekmer as skm  # kept for consistency
from snekmer.learn import generate_kmer_counts, match_kmer_counts_format


# ---------------------------------------------------------
# Core class
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
    annotation_files: List[str]
    input_data: str
    output_path: str
    annotation: List[pd.DataFrame]
    kmer_count_totals: pd.DataFrame
    seq_annot: Dict[str, str]
    kmer_list: List[str]
    seq_kmer_dict: Dict[str, List[int]]
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

    def construct_kmer_counts_dataframe(self) -> pd.DataFrame:
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
        return final_matrix_with_scores

    def filter_top_two_values(
        self, final_matrix_with_scores: pd.DataFrame
    ) -> pd.DataFrame:
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
        final_matrix_with_scores_write = pa.Table.from_pandas(
            final_matrix_with_scores
        )
        with gzip.open(self.output_path, "wb") as gzipped_file:
            csv.write_csv(final_matrix_with_scores_write, gzipped_file)

    def execute_all(
        self,
        save_apply_associations: bool = False,
        config: Optional[Mapping[str, Any]] = None,
    ) -> None:
        """
        Executes all the comparison steps in sequence.

        This includes:
            1. Loading inputs.
            2. Generating kmer counts.
            3. Adding known/unknown tags.
            4. Constructing a k-mer counts dataframe.
            5. Matching format with comparison data.
            6. Calculating cosine similarity.
            7. Optionally filtering to keep top two values.
            8. Writing results to output.

        Parameters
        ----------
        save_apply_associations : bool, optional
            If True, keep full similarity matrix.
            If False, keep only top two associations per sequence (set others to NaN).
            Ignored if `config` is provided.
        config : mapping, optional
            Full Snakemake config. If provided, uses
            config["learn_apply"]["save_apply_associations"]
            to decide whether to filter to top two.
        """
        # Decide behavior from config (Snakemake) or from the flag (demo/import use)
        if config is not None:
            try:
                save_apply_associations = bool(
                    config["learn_apply"]["save_apply_associations"]
                )
            except KeyError:
                # fall back to the default flag if the key is missing
                pass

        # 1–2: load inputs and build kmer counts
        self.generate_inputs()
        self.kmer_list, self.kmer_totals, self.seq_kmer_dict = generate_kmer_counts(
            self.input_data,
            self.kmer_list,
            self.kmer_totals,
            self.seq_kmer_dict,
            False,
        )

        # 3–4: tag and build kmer count table
        self.add_known_tag()
        kmer_counts = self.construct_kmer_counts_dataframe()

        # 5: match format with the compare_associations matrix
        kmer_counts, self.kmer_count_totals = match_kmer_counts_format(
            kmer_counts, self.kmer_count_totals
        )

        # 6: similarity
        final_matrix_with_scores = self.calculate_cosine_similarity(kmer_counts)

        # 7: optionally keep only top 2
        if not save_apply_associations:
            final_matrix_with_scores = self.filter_top_two_values(
                final_matrix_with_scores
            )

        # 8: write
        self.write_output(final_matrix_with_scores)


# ---------------------------------------------------------
# Snakemake entry point
# ---------------------------------------------------------

if "snakemake" in globals():
    cfg = snakemake.config
    analysis = EvaluateSequences(
        compare_associations=snakemake.input.compare_associations,
        annotation_files=list(snakemake.input.annotation),
        input_data=snakemake.input.data,
        output_path=snakemake.output.apply,
    )
    analysis.execute_all(config=cfg)
