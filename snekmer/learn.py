# snekmer/learn.py

# ---------------------------------------------------------
# Imports
# ---------------------------------------------------------

import itertools
import random
import sys

import numpy as np
import pandas as pd

from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import snekmer as skm

# ---------------------------------------------------------
# Classes / Functions
# ---------------------------------------------------------

def generate_kmer_counts(input_data, kmer_list: list, kmer_totals: list, seq_kmer_dict: dict, reverse: bool=False) -> None:
    """
    Generates k-mer counts for each sequence and aggregates totals.

    Args:
        input_data (str): Path to the input NPZ file containing sequence and k-mer list.
        kmer_list (list): List of k-mers to count.
        kmer_totals (list): List to accumulate total counts across all sequences.
        seq_kmer_dict (dict): Dictionary to populate with per-sequence k-mer count lists.
        reverse (bool): If True, reverse each sequence before counting.

    Returns:
        tuple: A tuple (kmer_list, kmer_totals, seq_kmer_dict) where:
            - kmer_list (list): The original list of k-mers.
            - kmer_totals (list): Updated k-mer total counts.
            - seq_kmer_dict (dict): Mapping from sequence ID to its k-mer count list.
    """
    kmer_list, df = skm.io.load_npz(input_data)
    kmer_list = kmer_list[0]
    seq_ids = df["sequence_id"]
    kmer_totals = [0] * len(kmer_list)
    kmer_len = len(kmer_list[0])

    for i, seq in enumerate(seq_ids):
        v = df["sequence"][i]
        if reverse:
            v = v[::-1]
        k_counts = {}
        items = [
            v[item : (item + kmer_len)]
            for item in range(0, (len((v)) - kmer_len + 1))
        ]
        for j in items:
            k_counts[j] = k_counts.get(j, 0) + 1
        store = [
            k_counts[item] if item in k_counts else 0
            for item in kmer_list
        ]
        for i, item in enumerate(kmer_list):
            if item in k_counts:
                kmer_totals[i] += k_counts[item]
        seq_kmer_dict[seq] = store
    return kmer_list, kmer_totals, seq_kmer_dict

def match_kmer_counts_format(kmer_counts: pd.DataFrame, kmer_count_totals: pd.DataFrame) -> Tuple[pd.DataFrame, 
                                                                                            pd.DataFrame]:
    """
    Matches the format of the provided kmer counts DataFrame to the format of the
    comparison data. Ensures columns align correctly.

    Args:
        kmer_counts (DataFrame): DataFrame of kmer counts to format.
        kmer_count_totals (DataFrame): Total Kmer Counts

    Returns:
        DataFrame: Formatted kmer counts DataFrame.
    """
    compare_check = len(str(kmer_counts.columns.values[10])) == len(str(kmer_count_totals.columns.values[10]))
    
    if compare_check:
        check = len(kmer_counts.columns.values)
        alphabet_initial = set(
            itertools.chain(
                *[list(x) for x in kmer_counts.columns.values[10:check]]
            )
        )
        alphabet_compare = set(
            itertools.chain(
                *[
                    list(x)
                    for x in kmer_count_totals.columns.values[10:check]
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
    kmer_count_totals.drop("Totals", axis=0, inplace=True)
    kmer_count_totals.drop("Kmer Count", axis=1, inplace=True)
    kmer_count_totals.drop("Sequence count", axis=1, inplace=True)

    column_order = list(
        set(kmer_counts.columns) | set(kmer_count_totals.columns)
    )
    kmer_counts = kmer_counts.reindex(columns=column_order, fill_value=0)
    kmer_count_totals = kmer_count_totals.reindex(
        columns=column_order, fill_value=0
    )

    return kmer_counts, kmer_count_totals


def apply_selection_method(
    seq_ann_scores: pd.DataFrame,
    selection_method: str,
    threshold_type: Optional[str],
    threshold_dict: dict,
    weight_top: Optional[str], 
    weight_distance: Optional[str]
) -> Tuple[List[Optional[str]], List[float], pd.DataFrame]:
    """
    Apply the chosen selection method to sequence annotation scores.

    Args:
        seq_ann_scores (pd.DataFrame): DataFrame of annotation scores.
        selection_method (str): The method to use ('top_hit', 'greatest_distance', etc.).
        threshold_type (Optional[str]): Which threshold to apply, or None.
        threshold_dict (dict): Mapping from annotation to threshold values.
        weight_top (Optional[str]): Weight for raw scores in combined method.
        weight_distance (Optional[str]): Weight for distances in combined method.

    Returns:
        Tuple[List[Optional[str]], List[float], pd.DataFrame]: predictions, deltas, top two scores.
    """
    if threshold_type == "None":
        threshold_type = None
    if selection_method == "top_hit" and threshold_type is None:
        return _select_top_hit(seq_ann_scores)
    if selection_method == "top_hit" and threshold_type is not None:
        return _select_top_hit_with_threshold(seq_ann_scores, threshold_dict)
    if selection_method == "greatest_distance":
        return _select_greatest_distance(seq_ann_scores, threshold_dict)
    if selection_method == "combined_distance":
        return _select_combined_distance(seq_ann_scores, threshold_dict, weight_top, weight_distance)
    else:
        raise ValueError(f"Invalid selection method: {selection_method}")

def _get_top_two(
    row: pd.Series,
    thresholds: Optional[np.ndarray] = None
) -> pd.Series:
    """
    Generic helper to extract the top two values (and their headers) from a row,
    optionally applying per-column thresholds before selection.
    """
    if thresholds is not None:
        row = row.where(row >= thresholds, np.nan)
    top2 = row.nlargest(2)
    return pd.Series({
        "key_value_one": top2.iloc[0] if len(top2) > 0 else np.nan,
        "key_value_one_header": top2.index[0] if len(top2) > 0 else np.nan,
        "key_value_two": top2.iloc[1] if len(top2) > 1 else np.nan,
        "key_value_two_header": top2.index[1] if len(top2) > 1 else np.nan,
    })

def _select_top_hit(
    seq_ann_scores: pd.DataFrame
) -> Tuple[List[Optional[str]], List[float], pd.DataFrame]:
    """
    Select the top hit based on raw cosine similarity scores.

    Args:
        seq_ann_scores (pd.DataFrame): DataFrame of cosine simialrity scores indexed by sequence.

    Returns:
        Tuple[List[Optional[str]], List[float], pd.DataFrame]:
            - predictions: list of top annotation per sequence,
            - deltas: list of differences between top two scores,
            - top_two: DataFrame of the top two scores for each sequence.
    """
    key_vals_df = seq_ann_scores.apply(
        lambda row: _get_top_two(row), axis=1
    )
    predictions = key_vals_df["key_value_one_header"].tolist()
    deltas = (key_vals_df["key_value_one"] - key_vals_df["key_value_two"]).tolist()
    top_two = key_vals_df[["key_value_one", "key_value_two"]]
    return predictions, deltas, top_two

def _select_top_hit_with_threshold(
    seq_ann_scores: pd.DataFrame,
    threshold_dict: dict
) -> Tuple[List[Optional[str]], List[float], pd.DataFrame]:
    """
    Select the top cosine similarity score (top hit) while applying thresholds.

    Args:
        seq_ann_scores (pd.DataFrame): DataFrame of annotation scores.
        threshold_dict (dict): Mapping from annotation to threshold value.

    Returns:
        Tuple[List[Optional[str]], List[float], pd.DataFrame]: predictions, deltas, and DataFrame of top two scores.
    """
    thresholds = seq_ann_scores.columns.to_series().map(
        threshold_dict
    ).to_numpy()
    key_vals_df = seq_ann_scores.apply(
        lambda row: _get_top_two(row, thresholds), axis=1
    )
    predictions = key_vals_df["key_value_one_header"].tolist()
    deltas = (key_vals_df["key_value_one"] - key_vals_df["key_value_two"]).tolist()
    top_two = key_vals_df[["key_value_one", "key_value_two"]]
    return predictions, deltas, top_two

def _select_greatest_distance(
    seq_ann_scores: pd.DataFrame,
    threshold_dict: dict
    ) -> Tuple[List[Optional[str]], List[float], pd.DataFrame]:
    """
    Select annotation based on greatest distance above threshold.

    Args:
        seq_ann_scores (pd.DataFrame): DataFrame of annotation scores.
        threshold_dict (dict): Mapping from annotation to threshold value.

    Returns:
        Tuple[List[Optional[str]], List[float], pd.DataFrame]:
            - predictions: list of annotations with greatest distance,
            - deltas: list of distances above threshold,
            - top_two: DataFrame of the two highest raw scores per sequence.
    """
    thresholds = seq_ann_scores.columns.to_series().map(threshold_dict)
    distances = seq_ann_scores.subtract(thresholds, axis=1)
    filtered_distances = distances.where(distances >= 0, np.nan)

    top_two_idx = np.argpartition(
        np.nan_to_num(filtered_distances.values, nan=-np.inf),
        -2,
        axis=1,
    )[:, -2:]

    sorted_idx = np.argsort(
        filtered_distances.values[np.arange(filtered_distances.shape[0])[:, None], top_two_idx],
        axis=1
    )[:, ::-1]

    idx_mat = top_two_idx[np.arange(filtered_distances.shape[0])[:, None], sorted_idx]
    top_two_scores = np.take_along_axis(seq_ann_scores.values, idx_mat, axis=1)
    top_two_distances = np.take_along_axis(filtered_distances.values, idx_mat, axis=1)
    top_two_headers = np.array(filtered_distances.columns)[idx_mat]

    flat_hdrs = top_two_headers.flatten()
    flat_th = thresholds[flat_hdrs].values
    top_thresh = flat_th.reshape(top_two_headers.shape)

    key_vals_df = pd.DataFrame({
        "key_value_one": top_two_scores[:, 0],
        "key_value_one_header": top_two_headers[:, 0],
        "key_value_one_distance": top_two_distances[:, 0],
        "key_value_one_threshold": top_thresh[:, 0],
        "key_value_two": top_two_scores[:, 1],
        "key_value_two_header": top_two_headers[:, 1],
        "key_value_two_distance": top_two_distances[:, 1],
        "key_value_two_threshold": top_thresh[:, 1],
    })

    deltas = key_vals_df["key_value_one_distance"].tolist()
    predictions = key_vals_df["key_value_one_header"].tolist()
    top_two = key_vals_df[["key_value_one", "key_value_two"]]
    return predictions, deltas, top_two

def _select_combined_distance(
    seq_ann_scores: pd.DataFrame,
    threshold_dict: dict,
    weight_top: str,
    weight_distance: str,
) -> Tuple[List[Optional[str]], List[float], pd.DataFrame]:
    """
    Select annotation using a combined metric of score and distance above threshold.

    Args:
        seq_ann_scores (pd.DataFrame): DataFrame of annotation scores.
        threshold_dict (dict): Mapping from annotation to threshold value.
        weight_top (float): Weight for raw scores.
        weight_distance (float): Weight for distance above threshold.

    Returns:
        Tuple[List[Optional[str]], List[float], pd.DataFrame]: predictions, deltas, and DataFrame of top two combined scores.
    """
    thresholds = seq_ann_scores.columns.to_series().map(threshold_dict)
    distances = seq_ann_scores - thresholds
    positive_distances = distances.where(distances >= 0, np.nan)
    combined_scores = (seq_ann_scores * weight_top) + (distances * weight_distance)
    combined_scores = combined_scores.where(positive_distances.notna(), np.nan)

    top_combined = combined_scores.max(axis=1)
    preds = combined_scores.idxmax(axis=1).where(~top_combined.isna(), None)

    deltas = []
    top_two_list = []

    # Method 1/2: Top hit with threshold
    filtered = seq_ann_scores.where(seq_ann_scores >= thresholds, np.nan)
    top1 = filtered.max(axis=1)
    fam1 = filtered.idxmax(axis=1)
    temp = filtered.apply(lambda row: row[row != row.max()], axis=1)
    second1 = temp.max(axis=1)

    # Method 3: Greatest distance from threshold
    pos_dist_3 = distances.where(distances >= 0, np.nan)
    top_dist_3 = pos_dist_3.max(axis=1)
    fam3 = pos_dist_3.idxmax(axis=1)

    for idx in range(len(preds)):
        pred = preds.iloc[idx]
        if pred is None:
            deltas.append(None)
            top_two_list.append([np.nan, np.nan])
            continue

        if pred == fam1.iloc[idx]:
            delta = top1.iloc[idx] - second1.iloc[idx]
            top_two_list.append([top1.iloc[idx], second1.iloc[idx]])
        elif pred == fam3.iloc[idx]:
            original = seq_ann_scores.iat[idx, seq_ann_scores.columns.get_loc(pred)]
            thr = thresholds[pred]
            delta = original - thr
            top_two_list.append([original, thr])
        else:
            row_comb = combined_scores.iloc[idx]
            second_comb = row_comb.drop(pred).max()
            delta = top_combined.iloc[idx] - second_comb
            top_two_list.append([top_combined.iloc[idx], second_comb])
        deltas.append(delta)

    top_two = pd.DataFrame(top_two_list, columns=["key_value_one", "key_value_two"])
    return preds.tolist(), deltas, top_two

def fragment(
    sequence: str,
    version: str,
    length: int,
    location: str,
    min_length: int
    ) -> List[str]:
    """
    Fragment a given sequence.

    Args:
        sequence (str): the sequence to fragment.
        version (str): fragmentation method, either "absolute" or "percent".
        length (int): length for fragmentation. If version is "percent", this is treated as a percentage.
        location (str): where the fragmentation happens - "start", "end", or "random".
        min_length (int): minimum length a fragment should be to be retained.

    Returns:
        List[str]: the retained fragment(s).
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
