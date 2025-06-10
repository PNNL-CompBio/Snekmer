# snekmer/learn.py

import numpy as np
import pandas as pd
from typing import Sequence, Dict, List, Mapping, Optional, Tuple
import itertools
import sys
import snekmer as skm



def generateKmerCounts(inputData, kmerList: list, kmerTotals: list, seqKmerdict: dict, reverse: bool=False) -> None:
    """
    Generates a dictionary of kmer counts for each sequence.

    Processes the input data to count the occurrence of each kmer in each sequence.
    """
    kmerList, df = skm.io.load_npz(inputData)
    kmerList = kmerList[0]
    seqIDs = df["sequence_id"]
    kmerTotals = [0] * len(kmerList)
    kmerLen = len(kmerList[0])

    for i, seq in enumerate(seqIDs):
        v = df["sequence"][i]
        if reverse:
            v = v[::-1]
        kCounts = {}
        items = [
            v[item : (item + kmerLen)]
            for item in range(0, (len((v)) - kmerLen + 1))
        ]
        for j in items:
            kCounts[j] = kCounts.get(j, 0) + 1
        store = [
            kCounts[item] if item in kCounts else 0
            for item in kmerList
        ]
        for i, item in enumerate(kmerList):
            if item in kCounts:
                kmerTotals[i] += kCounts[item]
        seqKmerdict[seq] = store
    return kmerList, kmerTotals, seqKmerdict

def matchKmerCountsFormat(kmerCounts: pd.DataFrame, kmerCountTotals: pd.DataFrame) -> Tuple[pd.DataFrame, 
                                                                                            pd.DataFrame]:
    """
    Matches the format of the provided kmer counts DataFrame to the format of the
    comparison data. Ensures columns align correctly.

    Args:
        kmerCounts (DataFrame): DataFrame of kmer counts to format.
        kmerCountTotals (DataFrame): Total Kmer Counts

    Returns:
        DataFrame: Formatted kmer counts DataFrame.
    """
    if len(str(kmerCounts.columns.values[10])) == len(
        str(kmerCountTotals.columns.values[10])
    ):
        compareCheck = True
    else:
        compareCheck = False

    if compareCheck:
        check = len(kmerCounts.columns.values)
        alphabetInitial = set(
            itertools.chain(
                *[list(x) for x in kmerCounts.columns.values[10:check]]
            )
        )
        alphabetCompare = set(
            itertools.chain(
                *[
                    list(x)
                    for x in kmerCountTotals.columns.values[10:check]
                ]
            )
        )
        if alphabetCompare != alphabetInitial:
            compareCheck = False

    if not compareCheck:
        print("Compare Check Failed. ")
        sys.exit()

    kmerCounts.drop("Totals", axis=0, inplace=True)
    kmerCounts.drop("Sequence count", axis=1, inplace=True)
    kmerCountTotals.drop("Totals", axis=0, inplace=True)
    kmerCountTotals.drop("Kmer Count", axis=1, inplace=True)
    kmerCountTotals.drop("Sequence count", axis=1, inplace=True)

    columnOrder = list(
        set(kmerCounts.columns) | set(kmerCountTotals.columns)
    )
    kmerCounts = kmerCounts.reindex(columns=columnOrder, fill_value=0)
    kmerCountTotals = kmerCountTotals.reindex(
        columns=columnOrder, fill_value=0
    )

    return kmerCounts, kmerCountTotals


def applySelectionMethod(
    seqAnnScores: pd.DataFrame,
    selectionMethod: str,
    thresholdType: Optional[str],
    thresholdDict: dict,
    weightTop: Optional[str], 
    weightDistance: Optional[str]
) -> Tuple[List[Optional[str]], List[float], pd.DataFrame]:
    """
    Apply the chosen selection method to sequence annotation scores.
    """
    if thresholdType == "None":
        thresholdType = None
    if selectionMethod == "top_hit" and thresholdType is None:
        return _select_top_hit(seqAnnScores)
    if selectionMethod == "top_hit" and thresholdType is not None:
        return _select_top_hit_with_threshold(seqAnnScores, thresholdDict)
    if selectionMethod == "greatest_distance":
        return _select_greatest_distance(seqAnnScores, thresholdDict)
    if selectionMethod == "combined_distance":
        return _select_combined_distance(seqAnnScores, thresholdDict, weightTop, weightDistance)
    else:
        raise ValueError(f"Invalid selection method: {selectionMethod}")

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
        "keyValueOne": top2.iloc[0] if len(top2) > 0 else np.nan,
        "keyValueOneHeader": top2.index[0] if len(top2) > 0 else np.nan,
        "keyValueTwo": top2.iloc[1] if len(top2) > 1 else np.nan,
        "keyValueTwoHeader": top2.index[1] if len(top2) > 1 else np.nan,
    })

def _select_top_hit(
    seqAnnScores: pd.DataFrame
) -> Tuple[List[Optional[str]], List[float], pd.DataFrame]:
    keyValsDf = seqAnnScores.apply(
        lambda row: _get_top_two(row), axis=1
    )
    predictions = keyValsDf["keyValueOneHeader"].tolist()
    deltas = (keyValsDf["keyValueOne"] - keyValsDf["keyValueTwo"]).tolist()
    topTwo = keyValsDf[["keyValueOne", "keyValueTwo"]]
    return predictions, deltas, topTwo

def _select_top_hit_with_threshold(
    seqAnnScores: pd.DataFrame,
    thresholdDict: dict
) -> Tuple[List[Optional[str]], List[float], pd.DataFrame]:
    thresholds = seqAnnScores.columns.to_series().map(
        thresholdDict
    ).to_numpy()
    keyValsDf = seqAnnScores.apply(
        lambda row: _get_top_two(row, thresholds), axis=1
    )
    predictions = keyValsDf["keyValueOneHeader"].tolist()
    deltas = (keyValsDf["keyValueOne"] - keyValsDf["keyValueTwo"]).tolist()
    topTwo = keyValsDf[["keyValueOne", "keyValueTwo"]]
    return predictions, deltas, topTwo

def _select_greatest_distance(
    seqAnnScores: pd.DataFrame,
    thresholdDict: dict
    ) -> Tuple[List[Optional[str]], List[float], pd.DataFrame]:
    
    thresholds = seqAnnScores.columns.to_series().map(thresholdDict)
    distances = seqAnnScores.subtract(thresholds, axis=1)
    filteredDistances = distances.where(distances >= 0, np.nan)

    topTwoIdx = np.argpartition(
        np.nan_to_num(filteredDistances.values, nan=-np.inf),
        -2,
        axis=1,
    )[:, -2:]

    sortedIdx = np.argsort(
        filteredDistances.values[np.arange(filteredDistances.shape[0])[:, None], topTwoIdx],
        axis=1
    )[:, ::-1]

    idx_mat = topTwoIdx[np.arange(filteredDistances.shape[0])[:, None], sortedIdx]
    topTwoScores = np.take_along_axis(seqAnnScores.values, idx_mat, axis=1)
    topTwoDistances = np.take_along_axis(filteredDistances.values, idx_mat, axis=1)
    topTwoHeaders = np.array(filteredDistances.columns)[idx_mat]

    flatHdrs = topTwoHeaders.flatten()
    flatTh = thresholds[flatHdrs].values
    topThresh = flatTh.reshape(topTwoHeaders.shape)

    keyValsDf = pd.DataFrame({
        "keyValueOne": topTwoScores[:, 0],
        "keyValueOneHeader": topTwoHeaders[:, 0],
        "keyValueOneDistance": topTwoDistances[:, 0],
        "keyValueOneThreshold": topThresh[:, 0],
        "keyValueTwo": topTwoScores[:, 1],
        "keyValueTwoHeader": topTwoHeaders[:, 1],
        "keyValueTwoDistance": topTwoDistances[:, 1],
        "keyValueTwoThreshold": topThresh[:, 1],
    })

    deltas = keyValsDf["keyValueOneDistance"].tolist()
    predictions = keyValsDf["keyValueOneHeader"].tolist()
    topTwo = keyValsDf[["keyValueOne", "keyValueTwo"]]
    return predictions, deltas, topTwo

def _select_combined_distance(
    seqAnnScores: pd.DataFrame,
    thresholdDict: dict,
    weightTop: str,
    weightDistance: str,
) -> Tuple[List[Optional[str]], List[float], pd.DataFrame]:
    thresholds = seqAnnScores.columns.to_series().map(thresholdDict)
    distances = seqAnnScores - thresholds
    positiveDistances = distances.where(distances >= 0, np.nan)
    combinedScores = (seqAnnScores * weightTop) + (distances * weightDistance)
    combinedScores = combinedScores.where(positiveDistances.notna(), np.nan)

    topCombined = combinedScores.max(axis=1)
    preds = combinedScores.idxmax(axis=1).where(~topCombined.isna(), None)

    deltas = []
    topTwoList = []

    # Method 1/2: Top hit with threshold
    filtered = seqAnnScores.where(seqAnnScores >= thresholds, np.nan)
    top1 = filtered.max(axis=1)
    fam1 = filtered.idxmax(axis=1)
    temp = filtered.apply(lambda row: row[row != row.max()], axis=1)
    second1 = temp.max(axis=1)

    # Method 3: Greatest distance from threshold
    posDist3 = distances.where(distances >= 0, np.nan)
    topDist3 = posDist3.max(axis=1)
    fam3 = posDist3.idxmax(axis=1)

    for idx in range(len(preds)):
        pred = preds.iloc[idx]
        if pred is None:
            deltas.append(None)
            topTwoList.append([np.nan, np.nan])
            continue

        if pred == fam1.iloc[idx]:
            delta = top1.iloc[idx] - second1.iloc[idx]
            topTwoList.append([top1.iloc[idx], second1.iloc[idx]])
        elif pred == fam3.iloc[idx]:
            original = seqAnnScores.iat[idx, seqAnnScores.columns.get_loc(pred)]
            thr = thresholds[pred]
            delta = original - thr
            topTwoList.append([original, thr])
        else:
            rowComb = combinedScores.iloc[idx]
            secondComb = rowComb.drop(pred).max()
            delta = topCombined.iloc[idx] - secondComb
            topTwoList.append([topCombined.iloc[idx], secondComb])
        deltas.append(delta)

    topTwo = pd.DataFrame(topTwoList, columns=["keyValueOne", "keyValueTwo"])
    return preds.tolist(), deltas, topTwo
