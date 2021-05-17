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
from kmerfeatures.model import KmerScaler


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
        Feature matrix in the form of rows (sequences), columns
        (kmers), where values are kmer counts for each protein.
    metric : str
        Metric used for pairwise distance computation (see
            sklearn.metrics.pairwise documentation).

    Returns
    -------
    numpy.ndarray
        Square matrix with the similarity scores of pairwise
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
    See sklearn documentation for more details.

    Parameters
    ----------
    feature_matrix : numpy.ndarray
        Feature matrix in the form of rows (sequences), columns
        (kmers), where values are kmer counts for each protein.
    method : type
        Description of parameter `method`.
    **kwargs : dict
        Keyword arguments for clustering class.

    Returns
    -------
    numpy.ndarray
        Array of clustered labels.

    """
    if method == "agglomerative":
        clusters = AgglomerativeClustering(n_clusters=n_clusters).fit_predict(feature_matrix)

    return clusters


def _apply_probability(kmer, label, compare_label):
    """Compute binary probability of a label for a kmer.

    Parameters
    ----------
    kmer : int
        Kmer count (assumes binary presence/absence (0, 1)).
    label : str
        Label value for a given kmer.
    compare_label : str
        Label value to assess a given kmer against.

    Returns
    -------
    list
        Binary probability (0, 1) of valid kmer with given label
        for all kmers in the kmer count array.
            e.g. [0, 1, 1, 0] for an input kmer array of length 4

    """
    return (kmer == 1) * 1 * (label == compare_label)


def _binarize(feature_matrix):
    """Convert feature counts into binary presence/absence.

    Parameters
    ----------
    feature_matrix : array (list, numpy.ndarray, pandas.Series)
        Feature matrix, where each row represents a kmer and each
        column represents a sequence

    Returns
    -------
    numpy.ndarray
        Description of returned object.

    """
    return (feature_matrix > 0) * 1.0


def feature_class_probabilities(feature_matrix, labels, kmers=None, processes=2):
    """Calculate probabilities for features being in a defined class.

    Note: only coded to work for the binary case (2 classes).

    Parameters
    ----------
    feature_matrix : array (list, numpy.ndarray, pandas.Series)
        Feature matrix, where each row represents a kmer and each
        column represents a sequence
        In other words, len(feature_matrix.T) must equal len(labels)
    labels : list or numpy.ndarray or pandas.Series
        Class labels describing feature matrix.
        Must have as many entries as the number of feature columns.
    kmers : list or None (default: None)
        Optional list of kmer identifiers that map to kmer columns.
        len(kmers) must equal len(feature_matrix).
    processes : int
        Number of processes (for multiprocessing)

    Returns
    -------
    numpy.ndarray (df=False) or pandas.DataFrame (df=True)
        Array or DataFrame containing scoring and probability
        results for each family.

    """
    # ensure labels are given as a numpy array
    if isinstance(labels, (np.ndarray, list, pd.Series)):
        labels = np.array(labels, dtype=object)
    else:
        raise TypeError("Labels must be list- or array-like.")

    # kmer labels
    if kmers is None:
        kmers = np.array([n for n in range(len(feature_matrix))])
    if len(kmers) != len(feature_matrix):
        raise ValueError("Kmer array shape is mismatched.")

    # check that labels are the same size as number of examples
    if len(feature_matrix.T) != len(labels):
        raise ValueError("Input shapes are mismatched.")

    # convert the feature count matrix into binary presence/absence
    feature_matrix = _binarize(feature_matrix)

    # iterate through every feature in the input matrix
    n_sequences = [
        len(np.hstack(np.where(labels == l))) for l in np.unique(labels)
        ]
    results = {l: {} for l in labels}
    for i, l in enumerate(np.unique(labels)):
        presence, probability = list(), list()
        matching_label = np.equal(labels, l)

        # for each row in feature matrix, compute each value's probability
        # no multiprocessing-- using numpy native parallelization
        for kmer in feature_matrix:
            matching_kmer = np.equal(kmer, 1)
            ones = np.ones(len(kmer))
            p = np.multiply(np.multiply(matching_kmer, matching_label), ones)
            presence.append(np.sum(p))
            probability.append(np.sum(p) / n_sequences[i])

        # if no kmers specified, save numerical range
        results[l]['kmer'] = kmers
        results[l]['count'] = np.asarray(presence, dtype=int)
        results[l]['probability'] = np.asarray(probability, dtype=float)

    # compute score based on (label, not label) assignment
    weight = 1 / (len(np.unique(labels)) - 1)
    for l in np.unique(labels):
        o = np.unique(labels)[np.unique(labels) != l][0]  # other labels
        results[l]['score'] = ((results[l]['probability'])
                               - (results[o]['probability'] * (weight)))

    results = pd.DataFrame(results).T.reset_index().rename(
        columns={'index': 'label'}
        )

    # reformat dataframe into long form
    results = results.set_index(['label']).apply(pd.Series.explode).reset_index()

    return results


# from Jason
# This code will apply the feature class probabilities derived in the previous functions
#   to a new set of examples as a weighted score. The idea here is to use a very simple
#   approach to classification.
def apply_feature_probabilities(
    feature_matrix, scores, scaler=False, **kwargs
):
    """Calculate class probability scores based on kmer vectors.

    Parameters
    ----------
    feature_matrix : array (list, numpy.ndarray, pandas.Series)
        Feature matrix of shape (n, m), where each row represents
        a kmer and each column represents a sequence.
        In other words, len(feature_matrix.T) must equal len(labels)
    scores : array (list, numpy.ndarray, pandas.Series)
        Array of shape (m,) of kmer probability scores.
        Scores must be in the same kmer order as in `feature_matrix`.
    scaler : bool (default: False)
        If True, performs scaling to reduce kmer features.
    **kwargs : dict
        Keyword arguments for KmerScaler().

    Returns
    -------
    numpy.ndarray
        Feature matrix weighted

    """
    # optionally apply scaling
    if scaler:
        scaler = KmerScaler(**kwargs)
        scaler.fit(scores)
        scores = scaler.transform(scores)
        feature_matrix = scaler.transform(feature_matrix)

    feature_matrix = _binarize(feature_matrix)
    score_matrix = scores * feature_matrix
    score_totals = np.array([np.sum(arr) for arr in score_matrix])

    return score_totals
