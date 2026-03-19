# ---------------------------------------------------------
# Imports
# ---------------------------------------------------------

import gzip
import json
import random
from typing import Any, Dict, List, Optional

import numpy as np
import pandas as pd

import snekmer as skm


# ---------------------------------------------------------
# Core functions
# ---------------------------------------------------------

def load_stats_from_csv(csv: str) -> Dict[str, Dict[str, Any]]:
    """
    Load existing family statistics from a csv checkpoint file.

    The csv must contain the following columns:
      - family (str): family identifier
      - count (int): number of observations
      - sum (float): sum of values
      - sumSqr (float): sum of squared values
      - min (float): minimum value
      - max (float): maximum value
      - percentileValues (str): JSON-encoded list of percentile thresholds

    Args:
        csv (str): Path to the csv file containing family statistics.

    Returns:
        Dict[str, Dict[str, Any]]:
    """
    csv = str(csv)
    if not csv:
        return {}

    df = pd.read_csv(csv)
    combined_stats: Dict[str, Dict[str, Any]] = {}

    families = df["family"].values
    counts = df["count"].values
    sums = df["sum"].values
    sum_sqrs = df["sumSqr"].values
    mins = df["min"].values
    maxs = df["max"].values
    pv_strings = df["percentileValues"].values

    for i in range(len(df)):
        combined_stats[families[i]] = {
            "count": int(counts[i]),
            "sum": float(sums[i]),
            "sumSqr": float(sum_sqrs[i]),
            "min": float(mins[i]),
            "max": float(maxs[i]),
            "percentileValues": json.loads(pv_strings[i]),
        }

    return combined_stats


def save_stats(
    combined_stats: Dict[str, Dict[str, Any]],
    csv: str,
) -> None:
    """
    Save combined family statistics to a csv checkpoint file.

    The output csv will contain the following columns:
      - family (str): the family identifier
      - count (int): number of observations
      - sum (float): sum of values
      - sumSqr (float): sum of squared values
      - min (float): minimum value
      - max (float): maximum value
      - percentileValues (str): JSON-encoded list of percentile thresholds

    Args:
        combined_stats (Dict[str, Dict[str, Any]]): Mapping from each family name to its statistics.
        csv (str): Path where the csv file will be written.
    """

    rows = []
    for family, stats in combined_stats.items():
        row = {
            "family": family,
            "count": stats["count"],
            "sum": stats["sum"],
            "sumSqr": stats["sumSqr"],
            "min": stats["min"],
            "max": stats["max"],
            "percentileValues": json.dumps(stats["percentileValues"]),
        }
        rows.append(row)

    df = pd.DataFrame(rows)
    df.to_csv(csv, index=False)


def _count_data_rows(filename: str) -> int:
    """Count data rows in a CSV file (gzipped or plain), excluding the header."""
    open_fn = gzip.open if filename.endswith(".gz") else open
    with open_fn(filename, "rt") as f:
        # Count lines and subtract 1 for header
        return sum(1 for _ in f) - 1


def collect_all_family_statistics(
    filenames: List[str],
    existing_stats: Optional[Dict[str, Dict[str, Any]]] = None,
    reservoir_size: int = 100000,
) -> Dict[str, Dict[str, Any]]:
    """
    Aggregate family-wise statistics from multiple CSV files using pre-determined
    random sampling for percentiles.

    Since all files are known upfront, this replaces reservoir sampling with a
    simple random sample: we count total rows across all files, pre-select which
    rows will be kept, then extract only those rows during processing. This is
    both faster and statistically cleaner than streaming reservoir sampling.

    For each family column, this function computes:
      - count: total number of values seen
      - sum: sum of all values
      - sumSqr: sum of squares of all values
      - min: minimum value
      - max: maximum value
      - percentileValues: a uniformly random sample of values (for approximate percentiles)

    Args:
        filenames (List[str]): Paths to the input CSV files (may be gzipped).
        existing_stats (Optional[Dict[str, Dict[str, Any]]]):
            Existing stats from a prior checkpoint. If None, starts fresh.
        reservoir_size (int): Maximum number of values to retain for percentile estimation.

    Returns:
        Dict[str, Dict[str, Any]]: Updated statistics per family.
    """
    if not filenames:
        return existing_stats if existing_stats else {}

    # ------------------------------------------------------------------
    # Phase 1: Count rows per file (fast — just counting newlines)
    # ------------------------------------------------------------------
    row_counts = [_count_data_rows(f) for f in filenames]
    total_new_rows = sum(row_counts)
    cumulative_rows = np.concatenate([[0], np.cumsum(row_counts)])

    # ------------------------------------------------------------------
    # Phase 2: Determine how many reservoir slots go to prior vs new data
    # ------------------------------------------------------------------
    prior_count = 0
    if existing_stats:
        first_fam = next(iter(existing_stats))
        prior_count = existing_stats[first_fam]["count"]

    grand_total = prior_count + total_new_rows

    if grand_total <= reservoir_size:
        # Everything fits — keep all new rows and all prior reservoir values
        n_from_new = total_new_rows
        n_from_prior = min(prior_count, reservoir_size) if existing_stats else 0
    else:
        # Proportional allocation: each item has equal probability of being selected
        n_from_new = int(round(reservoir_size * total_new_rows / grand_total))
        n_from_prior = reservoir_size - n_from_new

    # ------------------------------------------------------------------
    # Phase 3: Pre-select which new rows to keep (global indices)
    # ------------------------------------------------------------------
    if n_from_new >= total_new_rows:
        selected_new = np.arange(total_new_rows)
    else:
        selected_new = np.sort(
            np.random.choice(total_new_rows, n_from_new, replace=False)
        )

    # Map global new-data indices to per-file local indices
    file_selections = []
    for i in range(len(filenames)):
        start = cumulative_rows[i]
        end = cumulative_rows[i + 1]
        mask = (selected_new >= start) & (selected_new < end)
        local_indices = selected_new[mask] - start
        file_selections.append(local_indices)

    # ------------------------------------------------------------------
    # Phase 4: Process each file — bulk stats + extract pre-selected rows
    # ------------------------------------------------------------------
    family_names: Optional[List[str]] = None
    n_families: int = 0
    _sums: Optional[np.ndarray] = None
    _sum_sqrs: Optional[np.ndarray] = None
    _mins: Optional[np.ndarray] = None
    _maxs: Optional[np.ndarray] = None

    reservoir_parts: List[np.ndarray] = []

    for fi, filename in enumerate(filenames):
        df = pd.read_csv(filename, engine="c")
        current_families = df.columns[:-1].tolist()

        if _sums is None:
            family_names = current_families
            n_families = len(family_names)
        else:
            # Verify column order matches the first file
            if current_families != family_names:
                raise ValueError(
                    f"Column mismatch in {filename}: expected {len(family_names)} "
                    f"families in the same order as the first file, but got "
                    f"{len(current_families)} families. First difference at index "
                    f"{next(i for i, (a, b) in enumerate(zip(family_names, current_families)) if a != b) if len(family_names) == len(current_families) else 'length mismatch'}."
                )

        vals = df.iloc[:, :-1].to_numpy(dtype=np.float64)

        if _sums is None:
            _sums = np.zeros(n_families, dtype=np.float64)
            _sum_sqrs = np.zeros(n_families, dtype=np.float64)
            _mins = np.full(n_families, np.inf, dtype=np.float64)
            _maxs = np.full(n_families, -np.inf, dtype=np.float64)

            # Incorporate prior stats if present
            if existing_stats:
                for idx, fam in enumerate(family_names):
                    if fam in existing_stats:
                        s = existing_stats[fam]
                        _sums[idx] = s["sum"]
                        _sum_sqrs[idx] = s["sumSqr"]
                        _mins[idx] = s["min"]
                        _maxs[idx] = s["max"]

        # Bulk stats on ALL rows (4 numpy calls, no loops)
        _sums += vals.sum(axis=0)
        _sum_sqrs += np.einsum("ij,ij->j", vals, vals)
        np.minimum(_mins, vals.min(axis=0), out=_mins)
        np.maximum(_maxs, vals.max(axis=0), out=_maxs)

        # Extract only the pre-selected rows for the reservoir
        sel = file_selections[fi]
        if len(sel) > 0:
            reservoir_parts.append(vals[sel, :])

    # ------------------------------------------------------------------
    # Phase 5: Assemble final reservoir
    # ------------------------------------------------------------------
    if reservoir_parts:
        new_reservoir = np.vstack(reservoir_parts)  # (n_from_new, n_families)
    else:
        new_reservoir = np.empty((0, n_families), dtype=np.float64)

    # Merge with prior reservoir if checkpoint exists
    if existing_stats and n_from_prior > 0 and family_names:
        prior_res_len = len(existing_stats[family_names[0]]["percentileValues"])

        if prior_res_len > 0:
            if n_from_prior >= prior_res_len:
                prior_indices = np.arange(prior_res_len)
            else:
                prior_indices = np.sort(
                    np.random.choice(prior_res_len, n_from_prior, replace=False)
                )

            # Build prior reservoir as 2D array and subsample
            prior_reservoir = np.column_stack(
                [
                    np.array(existing_stats[fam]["percentileValues"])[prior_indices]
                    for fam in family_names
                ]
            )  # (n_from_prior, n_families)

            final_reservoir = np.vstack([prior_reservoir, new_reservoir])
        else:
            final_reservoir = new_reservoir
    else:
        final_reservoir = new_reservoir

    # Safety trim (shouldn't be needed, but guard against rounding)
    if final_reservoir.shape[0] > reservoir_size:
        final_reservoir = final_reservoir[:reservoir_size]

    # ------------------------------------------------------------------
    # Phase 6: Convert to dict format
    # ------------------------------------------------------------------
    total_count = prior_count + total_new_rows
    result: Dict[str, Dict[str, Any]] = {}
    for idx, fam in enumerate(family_names):
        result[fam] = {
            "count": int(total_count),
            "sum": float(_sums[idx]),
            "sumSqr": float(_sum_sqrs[idx]),
            "min": float(_mins[idx]),
            "max": float(_maxs[idx]),
            "percentileValues": final_reservoir[:, idx].tolist(),
        }

    return result


def generate_family_statistics(
    combined_stats: Dict[str, Dict[str, Any]]
) -> pd.DataFrame:
    """
    Generate a summary DataFrame of descriptive statistics for each family.

    Args:
        combined_stats (Dict[str, Dict[str, Any]]):
            Mapping from family names to their aggregated stats, where each stats dict contains:
            - count (int): number of observations
            - sum (float): sum of all values
            - sumSqr (float): sum of squared values
            - min (float): minimum observed value
            - max (float): maximum observed value
            - percentileValues (List[float]): sampled values used for percentile estimation

    Returns:
        pd.DataFrame:
            A table with one row per family and the following columns:
            - Family: family identifier
            - Mean: arithmetic mean of the values
            - Std_Dev: sample standard deviation
            - Min: minimum observed value
            - 10th_Percentile, 25th_Percentile, Median (50th), 75th_Percentile, 90th_Percentile
            - Max: maximum observed value
            - 1_Std_Dev_Above, 1_Std_Dev_Below: mean ± 1·std
            - 2_Std_Dev_Above, 2_Std_Dev_Below: mean ± 2·std
    """

    stats_data = {
        "family": [],
        "Mean": [],
        "Std_Dev": [],
        "Min": [],
        "10th_Percentile": [],
        "25th_Percentile": [],
        "Median": [],
        "75th_Percentile": [],
        "90th_Percentile": [],
        "Max": [],
        "1_Std_Dev_Above": [],
        "1_Std_Dev_Below": [],
        "2_Std_Dev_Above": [],
        "2_Std_Dev_Below": [],
    }

    for family, stats in combined_stats.items():
        n = stats["count"]
        cur_sum = stats["sum"]
        sum_sqr = stats["sumSqr"]
        mean = cur_sum / n if n > 0 else 0.0
        variance = (sum_sqr - (cur_sum ** 2) / n) / (n - 1) if n > 1 else 0.0
        standard_deviation = np.sqrt(variance)

        values = np.array(stats["percentileValues"])
        if len(values) > 0:
            percentiles = np.percentile(values, [10, 25, 50, 75, 90])
        else:
            # If no values, fill with NaN
            percentiles = [np.nan] * 5

        stats_data["family"].append(family)
        stats_data["Mean"].append(round(mean, 3))
        stats_data["Std_Dev"].append(round(standard_deviation, 3))
        stats_data["Min"].append(round(stats["min"], 3))
        stats_data["10th_Percentile"].append(round(percentiles[0], 3))
        stats_data["25th_Percentile"].append(round(percentiles[1], 3))
        stats_data["Median"].append(round(percentiles[2], 3))
        stats_data["75th_Percentile"].append(round(percentiles[3], 3))
        stats_data["90th_Percentile"].append(round(percentiles[4], 3))
        stats_data["Max"].append(round(stats["max"], 3))
        stats_data["1_Std_Dev_Above"].append(round(mean + standard_deviation, 3))
        stats_data["1_Std_Dev_Below"].append(round(mean - standard_deviation, 3))
        stats_data["2_Std_Dev_Above"].append(round(mean + 2 * standard_deviation, 3))
        stats_data["2_Std_Dev_Below"].append(round(mean - 2 * standard_deviation, 3))

    return pd.DataFrame(stats_data)


def execute_all(
    base_family_checkpoint: Optional[str],
    eval_apply_data: List[str],
    family_stats_output: str,
    checkpoint_output: str,
) -> None:
    """
    Executes the family statistics pipeline, essentially gathering values for family specific thresholds.

    Steps:
      1. Optionally load an existing family-level checkpoint CSV.
      2. Count rows across all files and pre-select which rows to sample.
      3. Process each file: compute bulk stats, extract pre-selected rows.
      4. Generate a summary DataFrame of descriptive statistics per family.
      5. Write the summary to `family_stats_output`.
      6. Save the updated combined_stats back to `checkpoint_output`.

    Args:
        base_family_checkpoint (Optional[str]): Path to an existing checkpoint CSV, or None.
        eval_apply_data (List[str]): List of paths to data files to incorporate.
        family_stats_output (str): Path where the summary statistics CSV will be written.
        checkpoint_output (str): Path where the updated checkpoint CSV will be written.
    """
    if base_family_checkpoint:
        existing_stats: Optional[Dict[str, Dict[str, Any]]] = load_stats_from_csv(
            base_family_checkpoint
        )
    else:
        existing_stats = None

    combined_stats = collect_all_family_statistics(
        filenames=eval_apply_data,
        existing_stats=existing_stats,
    )

    family_statistics_df = generate_family_statistics(combined_stats)
    family_statistics_df.to_csv(family_stats_output, index=False)
    save_stats(combined_stats, checkpoint_output)


# ---------------------------------------------------------
# Snakemake entry point
# ---------------------------------------------------------

if "snakemake" in globals():
    execute_all(
        base_family_checkpoint=snakemake.input.base_family_checkpoint,
        eval_apply_data=snakemake.input.eval_apply_data,
        family_stats_output=snakemake.output.family_stats,
        checkpoint_output=snakemake.output.checkpoint,
    )