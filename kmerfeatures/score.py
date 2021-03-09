"""score: Scoring and similarity analysis for extracted features.

author: @christinehc / @biodataganache
"""
# imports
from itertools import (product, repeat)
from multiprocessing import Pool
import numpy as np
import pandas as pd

from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics.pairwise import pairwise_distances


# functions
def to_feature_matrix(array):
    """Create properly shaped feature matrix for kmer scoring.

    Parameters
    ----------
    array : numpy.ndarray or list or array-like
        2-D array-like collection of kmer vectors.
        The assumed format is rows = sequences, cols = kmers.

    Returns
    -------
    numpy.ndarray of numpy.ndarrays
        2D array version of the 2D array-like input.

    """
    return np.array([np.array(a, dtype=int) for a in array])


def connection_matrix_from_features(feature_matrix, metric="jaccard"):
    """Calculate similarities based on features output from main().

    Parameters
    ----------
    feature_matrix : numpy.ndarray
        feature_matrix in the form of rows (sequences), columns (kmers),
        where values are the counts of the kmers for each protein.
    metric : str
        description.

    Returns
    -------
    numpy.ndarray
        a square matrix with the similarity scores of pairwise
        relationships between proteins.

    """
    if metric == "jaccard":
        sim = 1 - pairwise_distances(feature_matrix, metric="hamming")  # .T?
    else:
        sim = pairwise_distances(feature_matrix, metric=metric)
    return sim


def cluster_feature_matrix(feature_matrix,
                           method="agglomerative",
                           n_clusters=2):
    """Calculate clusters based on the feature matrix.

    Note: sklearn has a wide range of clustering options.

    Parameters
    ----------
    feature_matrix : pandas.DataFrame
        Description of parameter `feature_matrix`.
    method : type
        Description of parameter `method`.
    **kwargs : dict
        Keyword arguments for clustering class.

    Returns
    -------
    type
        Description of returned object.

    """
    # another place to have HPC since the clustering can be computationally
    #   intensive. AgglomerativeClustering may not be the best candidate
    #   for parallelization though - not sure

    if method == "agglomerative":
        clusters = AgglomerativeClustering(n_clusters=n_clusters).fit_predict(feature_matrix)

    return clusters


def _apply_probability(kmer, label, compare_label, processes=2):
    """Compute binary probability of a label for a kmer.

    Parameters
    ----------
    kmer : int
        Kmer count (assumes binary presence/absence (0, 1)).
    label : str
        Label value for a given kmer.
    compare_label : str
        Label value to assess a given kmer against.
    processes : int
        Number of processes for multiprocessing.

    Returns
    -------
    list
        Binary probability (0, 1) of valid kmer with given label
        for all kmers in the kmer count array.
            e.g. [0, 1, 1, 0] for an input kmer array of length 4

    """
    return (kmer == 1) * 1 * (label == compare_label)


def feature_class_probabilities(feature_matrix, labels, df=True, processes=2):
    """Calculate probabilities for features being in a defined class.

    Note: only coded to work for the binary case (2 classes).

    Parameters
    ----------
    feature_matrix : type
        Feature matrix, where each row represents a kmer and each
        column represents a sequence
        In other words, len(feature_matrix.T) must equal len(labels)
    labels : list or numpy.ndarray or pandas.Series
        Class labels describing feature matrix.
        Must have as many entries as the number of feature columns.
    df : bool
        If True, returns output as a pandas DataFrame;
        if False, returns output as a dictionary (default: True).
    processes : int
        Number of processes (for multiprocessing)

    Returns
    -------
    type
        Description of returned object.

    """
    # ensure labels are given as a numpy array
    if isinstance(labels, (np.ndarray, list, pd.Series)):
        labels = np.array(labels, dtype=object)
    else:
        raise TypeError("Labels must be list- or array-like.")

    # check that labels are the same size as number of examples
    if len(feature_matrix.T) != len(labels):
        raise ValueError("Input shapes are mismatched.")

    # convert the feature count matrix into binary presence/absence
    feature_matrix = (feature_matrix > 0) * 1.0

    # iterate through every feature in the input matrix
    results = {l: {'n_sequences': len(np.hstack(np.where(labels == l)))}
               for l in np.unique(labels)}
    for l in np.unique(labels):
        presence, probability = list(), list()
        matching_label = np.equal(labels, l)

        # for each row in feature matrix, compute each value's probability
        # attempt to use multiprocessing -- very slow
        # for i, kmer in enumerate(feature_matrix):
            # with Pool(processes) as pool:
            #     # kmer, label, compare_label, processes
            #     p = pool.starmap(_apply_probability, zip(
            #         kmer,
            #         repeat(labels[i]),
            #         repeat(l),
            #         repeat(processes)
            #         )
            #                      )

        # for each row in feature matrix, compute each value's probability
        # no multiprocessing-- using numpy native parallelization
        for kmer in feature_matrix:
            matching_kmer = np.equal(kmer, 1)
            ones = np.ones(len(kmer))
            p = np.multiply(np.multiply(matching_kmer, matching_label), ones)
            # p = np.array([(kmer[i] == 1) * 1 * (labels[i] == l)
            #               for i in range(len(kmer))])
            presence.append(p)
            probability.append(np.sum(p) / results[l]['n_sequences'])
        results[l]['presence'] = np.asarray(presence, dtype=int)
        results[l]['probability'] = np.asarray(probability, dtype=float)

    # compute score
    weight = 1 / (len(np.unique(labels)) - 1)
    for l in np.unique(labels):
        o = np.unique(labels)[np.unique(labels) != l][0]  # other labels
        results[l]['score'] = ((results[l]['probability'])
                               - (results[o]['probability'] * (weight)))

    if df:
        return pd.DataFrame(results).T.reset_index().rename(
            columns={'index': 'label'})

    return results
