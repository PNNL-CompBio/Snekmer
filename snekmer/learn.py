# snekmer/learn.py

# ---------------------------------------------------------
# Imports
# ---------------------------------------------------------

import itertools
import random
import sys
from collections import Counter

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

    # --- Optimization: pre-build kmer→index lookup for O(1) position access ---
    kmer_index = {km: idx for idx, km in enumerate(kmer_list)}
    num_kmers = len(kmer_list)

    for i, seq in enumerate(seq_ids):
        v = df["sequence"][i]
        if reverse:
            v = v[::-1]

        # Count kmers using Counter (marginally faster than manual .get loop)
        items = (
            v[pos : pos + kmer_len]
            for pos in range(len(v) - kmer_len + 1)
        )
        k_counts = Counter(items)

        # --- Optimization: build store array and update totals in one pass
        #     over observed kmers only (typically much smaller than kmer_list) ---
        store = [0] * num_kmers
        for kmer, count in k_counts.items():
            idx = kmer_index.get(kmer)
            if idx is not None:
                store[idx] = count
                kmer_totals[idx] += count

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


def _vectorized_top_two(
    values: np.ndarray,
    columns: np.ndarray,
    mask: Optional[np.ndarray] = None,
) -> pd.DataFrame:
    """
    Vectorized extraction of top-two values and their column headers from a 2D array.

    Parameters
    ----------
    values : np.ndarray
        2D array of scores (n_sequences x n_families).
    columns : np.ndarray
        1D array of column names corresponding to axis=1 of values.
    mask : np.ndarray, optional
        Boolean mask of same shape as values. Where False, values are excluded
        from consideration (treated as -inf for ranking purposes).

    Returns
    -------
    pd.DataFrame
        DataFrame with columns: key_value_one, key_value_one_header,
        key_value_two, key_value_two_header.
    """
    work = values.copy().astype(float)
    # NaN must be replaced with -inf so argpartition handles them correctly
    # (nlargest in the original code skips NaN; argpartition does not)
    nan_in_input = np.isnan(work)
    work[nan_in_input] = -np.inf
    if mask is not None:
        work[~mask] = -np.inf

    n_rows, n_cols = work.shape

    if n_cols >= 2:
        # argpartition is O(n) per row vs O(n log n) for full argsort
        top2_idx = np.argpartition(-work, 2, axis=1)[:, :2]
        # Sort the two candidates so index 0 is the largest
        top2_vals = np.take_along_axis(work, top2_idx, axis=1)
        order = np.argsort(-top2_vals, axis=1)
        top2_idx = np.take_along_axis(top2_idx, order, axis=1)
        top2_vals = np.take_along_axis(work, top2_idx, axis=1)
    else:
        # Edge case: 1 column
        top2_idx = np.zeros((n_rows, 2), dtype=int)
        top2_vals = np.full((n_rows, 2), -np.inf)
        top2_vals[:, 0] = work[:, 0]

    top2_headers = columns[top2_idx]

    # Replace -inf sentinels with NaN for output
    top2_vals[top2_vals == -np.inf] = np.nan

    # If mask was applied, entries where top value is NaN mean no valid hit
    result = pd.DataFrame({
        "key_value_one": top2_vals[:, 0],
        "key_value_one_header": top2_headers[:, 0],
        "key_value_two": top2_vals[:, 1],
        "key_value_two_header": top2_headers[:, 1],
    })

    # Where the top value is NaN, set header to NaN too
    nan_mask = np.isnan(top2_vals[:, 0])
    result.loc[nan_mask, "key_value_one_header"] = np.nan
    nan_mask2 = np.isnan(top2_vals[:, 1])
    result.loc[nan_mask2, "key_value_two_header"] = np.nan

    return result


def _select_top_hit(
    seq_ann_scores: pd.DataFrame
) -> Tuple[List[Optional[str]], List[float], pd.DataFrame]:
    """
    Select the top hit based on raw cosine similarity scores.

    Args:
        seq_ann_scores (pd.DataFrame): DataFrame of cosine similarity scores indexed by sequence.

    Returns:
        Tuple[List[Optional[str]], List[float], pd.DataFrame]:
            - predictions: list of top annotation per sequence,
            - deltas: list of differences between top two scores,
            - top_two: DataFrame of the top two scores for each sequence.
    """
    # --- Optimization: vectorized top-two extraction instead of row-wise apply ---
    key_vals_df = _vectorized_top_two(
        seq_ann_scores.values,
        np.array(seq_ann_scores.columns),
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
    # --- Optimization: vectorized top-two with threshold mask ---
    thresholds = seq_ann_scores.columns.to_series().map(
        threshold_dict
    ).to_numpy().astype(float)

    # Build boolean mask: True where score >= threshold
    mask = seq_ann_scores.values >= thresholds[np.newaxis, :]

    key_vals_df = _vectorized_top_two(
        seq_ann_scores.values,
        np.array(seq_ann_scores.columns),
        mask=mask,
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

    cols = np.array(combined_scores.columns)
    comb_vals = combined_scores.values
    scores_vals = seq_ann_scores.values
    thresh_vals = thresholds.values.astype(float)
    dist_vals = distances.values

    n_rows = comb_vals.shape[0]

    # --- Optimization: vectorized prediction via argmax on nan-safe array ---
    comb_work = np.where(np.isnan(comb_vals), -np.inf, comb_vals)
    top_comb_idx = np.argmax(comb_work, axis=1)
    top_comb_val = comb_work[np.arange(n_rows), top_comb_idx]

    # Where all values were NaN, prediction is None
    all_nan = np.all(np.isnan(comb_vals), axis=1)
    preds_arr = cols[top_comb_idx]

    # Method 1: Top hit with threshold (filtered scores)
    filtered_vals = np.where(scores_vals >= thresh_vals[np.newaxis, :], scores_vals, -np.inf)
    top1_idx = np.argmax(filtered_vals, axis=1)
    top1_val = filtered_vals[np.arange(n_rows), top1_idx]
    fam1 = cols[top1_idx]

    # Second-best filtered score: mask out ALL values tied at the max
    # (original does row[row != row.max()].max() which excludes all ties)
    filtered_second = filtered_vals.copy()
    tied_at_max = (filtered_vals == top1_val[:, np.newaxis])
    filtered_second[tied_at_max] = -np.inf
    second1_val = np.max(filtered_second, axis=1)

    # Method 3: Greatest positive distance
    pos_dist_vals = np.where(dist_vals >= 0, dist_vals, -np.inf)
    top_dist_idx = np.argmax(pos_dist_vals, axis=1)
    fam3 = cols[top_dist_idx]

    # Compute deltas and top_two vectorized by case
    is_fam1 = (preds_arr == fam1)
    is_fam3 = (preds_arr == fam3) & ~is_fam1

    # Case 1: pred matches top filtered hit
    delta_case1 = top1_val - second1_val
    tv1_one = top1_val
    tv1_two = second1_val

    # Case 3: pred matches greatest distance family
    # Use dict lookup instead of searchsorted (columns are not sorted)
    col_to_idx = {c: i for i, c in enumerate(cols)}
    pred_col_idx = np.array([
        col_to_idx.get(p, 0) for p in preds_arr
    ])  # fallback index for None preds doesn't matter, will be masked
    original_scores = scores_vals[np.arange(n_rows), pred_col_idx]
    pred_thresholds = thresh_vals[pred_col_idx]
    delta_case3 = original_scores - pred_thresholds
    tv3_one = original_scores
    tv3_two = pred_thresholds

    # Case else: combined score difference
    comb_second = comb_work.copy()
    comb_second[np.arange(n_rows), top_comb_idx] = -np.inf
    second_comb_val = np.max(comb_second, axis=1)
    delta_else = top_comb_val - second_comb_val
    tve_one = top_comb_val
    tve_two = second_comb_val

    # Assemble results
    deltas_arr = np.where(is_fam1, delta_case1,
                 np.where(is_fam3, delta_case3, delta_else))
    tv_one = np.where(is_fam1, tv1_one,
             np.where(is_fam3, tv3_one, tve_one))
    tv_two = np.where(is_fam1, tv1_two,
             np.where(is_fam3, tv3_two, tve_two))

    # Handle None predictions
    preds_list = [None if all_nan[i] else preds_arr[i] for i in range(n_rows)]
    deltas_list = [None if all_nan[i] else float(deltas_arr[i]) for i in range(n_rows)]
    tv_one = np.where(all_nan, np.nan, tv_one)
    tv_two = np.where(all_nan, np.nan, tv_two)

    # Replace -inf sentinels with NaN
    tv_one = np.where(tv_one == -np.inf, np.nan, tv_one)
    tv_two = np.where(tv_two == -np.inf, np.nan, tv_two)

    top_two = pd.DataFrame({"key_value_one": tv_one, "key_value_two": tv_two})
    return preds_list, deltas_list, top_two

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