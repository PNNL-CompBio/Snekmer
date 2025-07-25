# ---------------------------------------------------------
# Imports
# ---------------------------------------------------------

import json
import random

import numpy as np
import pandas as pd

from typing import Any, Dict, List, Optional

import snekmer as skm

# ---------------------------------------------------------
# Files and Parameters
# ---------------------------------------------------------

config = snakemake.config

# ---------------------------------------------------------
# Run script
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
    combined_stats: dict[str, dict] = {}

    for _, row in df.iterrows():
        family = row["family"]
        combined_stats[family] = {
            "count": int(row["count"]),
            "sum": float(row["sum"]),
            "sumSqr": float(row["sumSqr"]),
            "min": float(row["min"]),
            "max": float(row["max"]),
            "percentileValues": json.loads(row["percentileValues"]),
        }

    return combined_stats

def save_stats(
    combined_stats: Dict[str, Dict[str, Any]],
    csv: str
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
            "percentileValues": json.dumps(
                stats["percentileValues"]
            ),
        }
        rows.append(row)

    df = pd.DataFrame(rows)
    df.to_csv(csv, index=False)


def collect_family_statistics(filename: str,
    existing_stats: Optional[Dict[str, Dict[str, Any]]] = None
    ) -> Dict[str, Dict[str, Any]]:
    """
    Aggregate family-wise statistics from a CSV file using reservoir sampling for percentiles.

    Reads the input CSV in chunks, updating for each family column:
      - count: total number of values seen
      - sum: sum of all values
      - sumSqr: sum of squares of all values
      - min: minimum observed value
      - max: maximum observed value
      - percentileValues: a reservoir‐sampled list of values (for approximate percentiles)

    Args:
        filename (str): Path to the input CSV file.
        existing_stats (Optional[Dict[str, Dict[str, Any]]]): 
            Existing stats to update. If None, a fresh stats dict is created.
    """
    
    chunk_size = 10000 
    reservoir_size = 100000

    if existing_stats is None:
        existing_stats = {}

    for chunk in pd.read_csv(filename, chunksize=chunk_size, engine="c"):
        families = chunk.columns[
            :-1
        ]  # Exclude the last column if it's the sequence name

        for family in families:
            values = chunk[family].dropna().astype(float).values
            if family not in existing_stats:
                existing_stats[family] = {
                    "count": 0,
                    "sum": 0.0,
                    "sumSqr": 0.0,
                    "min": np.inf,
                    "max": -np.inf,
                    "percentileValues": [],
                }

            stats = existing_stats[family]

            n = len(values)
            if n == 0:
                continue

                # Update count and sum statistics
            stats["sum"] += values.sum()
            stats["sumSqr"] += np.dot(values, values)
            stats["min"] = min(stats["min"], values.min())
            stats["max"] = max(stats["max"], values.max())

            # Reservoir sampling for percentiles
            for value in values:
                stats["count"] += 1  # Update total count
                total_seen = stats["count"]

                if len(stats["percentileValues"]) < reservoir_size:
                    # Fill the reservoir until it reaches the desired size
                    stats["percentileValues"].append(value)
                else:
                    # Replace elements with decreasing probability
                    j = random.randint(0, total_seen - 1)
                    if j < reservoir_size:
                        stats["percentileValues"][j] = value

        del chunk  # Free memory

    return existing_stats


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
            - Std Dev: sample standard deviation
            - Min: minimum observed value
            - 10th Percentile, 25th Percentile, Median (50th), 75th Percentile, 90th Percentile
            - Max: maximum observed value
            - 1 Std Dev Above, 1 Std Dev Below: mean ± 1·std
            - 2 Std Dev Above, 2 Std Dev Below: mean ± 2·std
    """
    
    stats_data = {
        "family": [],
        "Mean": [],
        "Std Dev": [],
        "Min": [],
        "10th Percentile": [],
        "25th Percentile": [],
        "Median": [],
        "75th Percentile": [],
        "90th Percentile": [],
        "Max": [],
        "1 Std Dev Above": [],
        "1 Std Dev Below": [],
        "2 Std Dev Above": [],
        "2 Std Dev Below": [],
    }

    for family, stats in combined_stats.items():
        n = stats["count"]
        cur_sum = stats["sum"]
        sum_sqr = stats["sumSqr"]
        mean = cur_sum / n if n > 0 else 0.0
        variance = (sum_sqr - (cur_sum**2) / n) / (n - 1) if n > 1 else 0.0
        standard_deviation = np.sqrt(variance)

        values = np.array(stats["percentileValues"])
        if len(values) > 0:
            percentiles = np.percentile(
                values, [10, 25, 50, 75, 90]
            )
        else:
            # If no values, fill with NaN
            percentiles = [np.nan] * 11

        stats_data["family"].append(family)
        stats_data["Mean"].append(round(mean, 3))
        stats_data["Std Dev"].append(round(standard_deviation, 3))
        stats_data["Min"].append(round(stats["min"], 3))
        stats_data["10th Percentile"].append(round(percentiles[0], 3))
        stats_data["25th Percentile"].append(round(percentiles[1], 3))
        stats_data["Median"].append(round(percentiles[2], 3))
        stats_data["75th Percentile"].append(round(percentiles[3], 3))
        stats_data["90th Percentile"].append(round(percentiles[4], 3))
        stats_data["Max"].append(round(stats["max"], 3))
        stats_data["1 Std Dev Above"].append(round(mean + standard_deviation, 3))
        stats_data["1 Std Dev Below"].append(round(mean - standard_deviation, 3))
        stats_data["2 Std Dev Above"].append(round(mean + 2 * standard_deviation, 3))
        stats_data["2 Std Dev Below"].append(round(mean - 2 * standard_deviation, 3))

    return pd.DataFrame(stats_data)

def execute_all(
    base_family_checkpoint: Optional[str],
    eval_apply_data: List[str],
    family_stats_output: str,
    checkpoint_output: str
) -> None:
    """
    Executes the family statistics pipeline, essentially gathering values for family specific thresholds.

    Steps:
      1. Optionally load an existing family‐level checkpoint CSV.
      2. Iterate over each evaluation/apply data file and update statistics via reservoir sampling.
      3. Generate a summary DataFrame of descriptive statistics per family.
      4. Write the summary to `family_stats_output`.
      5. Save the updated combined_stats back to `checkpoint_output`.

    Args:
        base_family_checkpoint (Optional[str]): Path to an existing checkpoint CSV, or None.
        eval_apply_data (List[str]): List of paths to data files to incorporate.
        family_stats_output (str): Path where the summary statistics CSV will be written.
        checkpoint_output (str): Path where the updated checkpoint CSV will be written.

    Returns:
        None
    """
    if base_family_checkpoint:
        combined_stats: Dict[str, Dict[str, Any]] = load_stats_from_csv(base_family_checkpoint)
    else:
        combined_stats = {}
    for filename in eval_apply_data:
        combined_stats = collect_family_statistics(
            filename,
            existing_stats=combined_stats
        )

    family_statistics_df = generate_family_statistics(combined_stats)
    family_statistics_df.to_csv(family_stats_output, index=False)
    save_stats(combined_stats, checkpoint_output)


execute_all(snakemake.input.base_family_checkpoint,
           snakemake.input.eval_apply_data,
           snakemake.output.family_stats,
           snakemake.output.checkpoint)
