# ---------------------------------------------------------
# Imports
# ---------------------------------------------------------

import pickle
from datetime import datetime
from os import makedirs
from os.path import exists, join

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import snekmer as skm
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split, StratifiedKFold
from sklearn.preprocessing import LabelEncoder

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.csv as csv
import seaborn as sns
import sklearn
from Bio import SeqIO
from scipy.interpolate import interp1d
from scipy.stats import rankdata
import random
import snekmer as skm
from scipy.ndimage import gaussian_filter1d
import sklearn.metrics.pairwise
import itertools
import os
import shutil
import sys
import json
# ---------------------------------------------------------
# Files and Parameters
# ---------------------------------------------------------

config = snakemake.config

# change matplotlib backend to non-interactive
plt.switch_backend("Agg")

# ---------------------------------------------------------
# Run script
# ---------------------------------------------------------


def load_existing_stats_from_csv(csv_file):
    """Load existing family statistics from a CSV checkpoint file.

    The CSV must contain the columns:
    Family, count, sum, sumSqr, min, max, values_for_percentiles
    (values_for_percentiles is a JSON-encoded list)
    """
    csv_file = str(csv_file)
    print(f"csv_file:{csv_file}")
    print(f"csv_file type:{type(csv_file)}")

    if not csv_file:
        return {}

    df = pd.read_csv(csv_file)
    combinedStats: dict[str, dict] = {}

    for _, row in df.iterrows():
        family = row["Family"]
        combinedStats[family] = {
            "count": int(row["count"]),
            "sum": float(row["sum"]),
            "sumSqr": float(row["sumSqr"]),
            "min": float(row["min"]),
            "max": float(row["max"]),
            "values_for_percentiles": json.loads(row["values_for_percentiles"]),
        }

    return combinedStats


def save_stats_to_csv(combinedStats, csv_file):

    rows = []
    for family, stats in combinedStats.items():
        row = {
            "Family": family,
            "count": stats["count"],
            "sum": stats["sum"],
            "sumSqr": stats["sumSqr"],
            "min": stats["min"],
            "max": stats["max"],
            # Encode values_for_percentiles as JSON string
            "values_for_percentiles": json.dumps(
                stats["values_for_percentiles"]
            ),
        }
        rows.append(row)

    df = pd.DataFrame(rows)
    df.to_csv(csv_file, index=False)


def collectFamilyStatistics(filename, existingStats=None):
    chunk_size = 10000  # Adjust based on available memory
    reservoir_size = 100000  # Size of the reservoir for percentiles

    if existingStats is None:
        existingStats = {}

    for chunk in pd.read_csv(filename, chunksize=chunk_size, engine="c"):
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
                    "values_for_percentiles": [],
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

                if len(stats["values_for_percentiles"]) < reservoir_size:
                    # Fill the reservoir until it reaches the desired size
                    stats["values_for_percentiles"].append(value)
                else:
                    # Replace elements with decreasing probability
                    j = random.randint(0, total_seen - 1)
                    if j < reservoir_size:
                        stats["values_for_percentiles"][j] = value

        del chunk  # Free memory

    return existingStats


def generateFamilyStatistics(combinedStats):
    statsData = {
        "Family": [],
        "Mean": [],
        "Std Dev": [],
        "Min": [],
        "10th Percentile": [],
        "20th Percentile": [],
        "25th Percentile": [],
        "30th Percentile": [],
        "40th Percentile": [],
        "Median": [],
        "60th Percentile": [],
        "70th Percentile": [],
        "75th Percentile": [],
        "80th Percentile": [],
        "90th Percentile": [],
        "Max": [],
        "1 Std Dev Above": [],
        "1 Std Dev Below": [],
        "2 Std Dev Above": [],
        "2 Std Dev Below": [],
    }

    for family, stats in combinedStats.items():
        n = stats["count"]
        sum_ = stats["sum"]
        sumSqr = stats["sumSqr"]
        mean = sum_ / n if n > 0 else 0.0
        variance = (sumSqr - (sum_**2) / n) / (n - 1) if n > 1 else 0.0
        std_dev = np.sqrt(variance)

        values = np.array(stats["values_for_percentiles"])
        if len(values) > 0:
            percentiles = np.percentile(
                values, [10, 20, 25, 30, 40, 50, 60, 70, 75, 80, 90]
            )
        else:
            # If no values, fill with NaN
            percentiles = [np.nan] * 11

        statsData["Family"].append(family)
        statsData["Mean"].append(round(mean, 3))
        statsData["Std Dev"].append(round(std_dev, 3))
        statsData["Min"].append(round(stats["min"], 3))
        statsData["10th Percentile"].append(round(percentiles[0], 3))
        statsData["20th Percentile"].append(round(percentiles[1], 3))
        statsData["25th Percentile"].append(round(percentiles[2], 3))
        statsData["30th Percentile"].append(round(percentiles[3], 3))
        statsData["40th Percentile"].append(round(percentiles[4], 3))
        statsData["Median"].append(round(percentiles[5], 3))
        statsData["60th Percentile"].append(round(percentiles[6], 3))
        statsData["70th Percentile"].append(round(percentiles[7], 3))
        statsData["75th Percentile"].append(round(percentiles[8], 3))
        statsData["80th Percentile"].append(round(percentiles[9], 3))
        statsData["90th Percentile"].append(round(percentiles[10], 3))
        statsData["Max"].append(round(stats["max"], 3))
        statsData["1 Std Dev Above"].append(round(mean + std_dev, 3))
        statsData["1 Std Dev Below"].append(round(mean - std_dev, 3))
        statsData["2 Std Dev Above"].append(round(mean + 2 * std_dev, 3))
        statsData["2 Std Dev Below"].append(round(mean - 2 * std_dev, 3))

    return pd.DataFrame(statsData)

    # Load existing stats from checkpoint CSV if exists


if snakemake.input.baseFamilyCheckpoint:
    print(f"input.baseFamilyCheckpoint is: {snakemake.input.baseFamilyCheckpoint}")
    combinedStats = load_existing_stats_from_csv(snakemake.input.baseFamilyCheckpoint)
else:
    combinedStats = None

    # Update combinedStats with new data
for filename in snakemake.input.evalApplyData:
    combinedStats = collectFamilyStatistics(
        filename, existingStats=combinedStats
    )

    # Generate updated family statistics
familyStatisticsDf = generateFamilyStatistics(combinedStats)
familyStatisticsDf.to_csv(snakemake.output.familyStats, index=False)

# Save updated combinedStats to checkpoint CSV
save_stats_to_csv(combinedStats, snakemake.output.checkpoint)