# ---------------------------------------------------------
# Imports
# ---------------------------------------------------------

import json
import random
import numpy as np
from typing import Dict, Any, List, Optional
import pandas as pd
import snekmer as skm

# ---------------------------------------------------------
# Files and Parameters
# ---------------------------------------------------------

config = snakemake.config

# ---------------------------------------------------------
# Run script
# ---------------------------------------------------------


def loadStatsFromCSV(csv: str) -> Dict[str, Dict[str, Any]]:    
    """
    Load existing family statistics from a csv checkpoint file.

    The csv must contain the following columns:
      - Family (str): family identifier
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
    combinedStats: dict[str, dict] = {}

    for _, row in df.iterrows():
        family = row["Family"]
        combinedStats[family] = {
            "count": int(row["count"]),
            "sum": float(row["sum"]),
            "sumSqr": float(row["sumSqr"]),
            "min": float(row["min"]),
            "max": float(row["max"]),
            "percentileValues": json.loads(row["percentileValues"]),
        }

    return combinedStats

def saveStats(
    combinedStats: Dict[str, Dict[str, Any]],
    csv: str
    ) -> None:
    """
    Save combined family statistics to a csv checkpoint file.

    The output csv will contain the following columns:
      - Family (str): the family identifier
      - count (int): number of observations
      - sum (float): sum of values
      - sumSqr (float): sum of squared values
      - min (float): minimum value
      - max (float): maximum value
      - percentileValues (str): JSON-encoded list of percentile thresholds

    Args:
        combinedStats (Dict[str, Dict[str, Any]]): Mapping from each family name to its statistics.
        csv (str): Path where the csv file will be written.
    """

    rows = []
    for family, stats in combinedStats.items():
        row = {
            "Family": family,
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


def collectFamilyStatistics(filename: str,
    existingStats: Optional[Dict[str, Dict[str, Any]]] = None
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
        existingStats (Optional[Dict[str, Dict[str, Any]]]): 
            Existing stats to update. If None, a fresh stats dict is created.
    """
    
    chunkSize = 10000 
    reservoirSize = 100000

    if existingStats is None:
        existingStats = {}

    for chunk in pd.read_csv(filename, chunksize=chunkSize, engine="c"):
        families = chunk.columns[
            :-1
        ]  # Exclude the last column if it's the sequence name

        for family in families:
            values = chunk[family].dropna().astype(float).values
            if family not in existingStats:
                existingStats[family] = {
                    "count": 0,
                    "sum": 0.0,
                    "sumSqr": 0.0,
                    "min": np.inf,
                    "max": -np.inf,
                    "percentileValues": [],
                }

            stats = existingStats[family]

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

                if len(stats["percentileValues"]) < reservoirSize:
                    # Fill the reservoir until it reaches the desired size
                    stats["percentileValues"].append(value)
                else:
                    # Replace elements with decreasing probability
                    j = random.randint(0, total_seen - 1)
                    if j < reservoirSize:
                        stats["percentileValues"][j] = value

        del chunk  # Free memory

    return existingStats


def generateFamilyStatistics(
        combinedStats: Dict[str, Dict[str, Any]]
    ) -> pd.DataFrame:
    """
    Generate a summary DataFrame of descriptive statistics for each family.

    Args:
        combinedStats (Dict[str, Dict[str, Any]]):
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
    
    statsData = {
        "Family": [],
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

    for family, stats in combinedStats.items():
        n = stats["count"]
        curSum = stats["sum"]
        sumSqr = stats["sumSqr"]
        mean = curSum / n if n > 0 else 0.0
        variance = (sumSqr - (curSum**2) / n) / (n - 1) if n > 1 else 0.0
        standardDeviation = np.sqrt(variance)

        values = np.array(stats["percentileValues"])
        if len(values) > 0:
            percentiles = np.percentile(
                values, [10, 25, 50, 75, 90]
            )
        else:
            # If no values, fill with NaN
            percentiles = [np.nan] * 11

        statsData["Family"].append(family)
        statsData["Mean"].append(round(mean, 3))
        statsData["Std Dev"].append(round(standardDeviation, 3))
        statsData["Min"].append(round(stats["min"], 3))
        statsData["10th Percentile"].append(round(percentiles[0], 3))
        statsData["25th Percentile"].append(round(percentiles[1], 3))
        statsData["Median"].append(round(percentiles[2], 3))
        statsData["75th Percentile"].append(round(percentiles[3], 3))
        statsData["90th Percentile"].append(round(percentiles[4], 3))
        statsData["Max"].append(round(stats["max"], 3))
        statsData["1 Std Dev Above"].append(round(mean + standardDeviation, 3))
        statsData["1 Std Dev Below"].append(round(mean - standardDeviation, 3))
        statsData["2 Std Dev Above"].append(round(mean + 2 * standardDeviation, 3))
        statsData["2 Std Dev Below"].append(round(mean - 2 * standardDeviation, 3))

    return pd.DataFrame(statsData)


# # Run all

# if snakemake.input.baseFamilyCheckpoint:
#     print(f"input.baseFamilyCheckpoint is: {snakemake.input.baseFamilyCheckpoint}")
#     combinedStats = loadStatsFromCSV(snakemake.input.baseFamilyCheckpoint)
# else:
#     combinedStats = None

#     # Update combinedStats with new data
# for filename in snakemake.input.evalApplyData:
#     combinedStats = collectFamilyStatistics(
#         filename, existingStats=combinedStats
#     )

#     # Generate updated family statistics
# familyStatisticsDf = generateFamilyStatistics(combinedStats)
# familyStatisticsDf.to_csv(snakemake.output.familyStats, index=False)

# # Save updated combinedStats to checkpoint csv
# saveStats(combinedStats, snakemake.output.checkpoint)


def executeAll(
    baseFamilyCheckpoint: Optional[str],
    evalApplyData: List[str],
    familyStatsOutput: str,
    checkpointOutput: str
) -> None:
    """
    Executes the family statistics pipeline, essentially gathering values for family specific thresholds.

    Steps:
      1. Optionally load an existing family‐level checkpoint CSV.
      2. Iterate over each evaluation/apply data file and update statistics via reservoir sampling.
      3. Generate a summary DataFrame of descriptive statistics per family.
      4. Write the summary to `familyStatsOutput`.
      5. Save the updated combinedStats back to `checkpointOutput`.

    Args:
        baseFamilyCheckpoint (Optional[str]): Path to an existing checkpoint CSV, or None.
        evalApplyData (List[str]): List of paths to data files to incorporate.
        familyStatsOutput (str): Path where the summary statistics CSV will be written.
        checkpointOutput (str): Path where the updated checkpoint CSV will be written.

    Returns:
        None
    """
    # 1. Load or initialize combinedStats
    if baseFamilyCheckpoint:
        combinedStats: Dict[str, Dict[str, Any]] = loadStatsFromCSV(baseFamilyCheckpoint)
    else:
        combinedStats = {}

    # 2. Update stats with each data file
    for filename in evalApplyData:
        combinedStats = collectFamilyStatistics(
            filename,
            existingStats=combinedStats
        )

    # 3. Generate descriptive statistics DataFrame
    familyStatisticsDf = generateFamilyStatistics(combinedStats)

    # 4. Write summary statistics
    familyStatisticsDf.to_csv(familyStatsOutput, index=False)

    # 5. Save updated checkpoint
    saveStats(combinedStats, checkpointOutput)


executeAll(snakemake.input.baseFamilyCheckpoint,
           snakemake.input.evalApplyData,
           snakemake.output.familyStats,
           snakemake.output.checkpoint)