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
from .model import KmerScoreScaler


# classes
class KmerScorer:
    """Score kmer vectors based on summed probability scores.

    Attributes
    ----------
    probabilities : numpy.ndarray
        Array of shape (n,) containing probabilities for n kmers.
    kmers : numpy.ndarray
        Array of shape (n, ) containing labels for n kmers.

    """
    def __init__(self):
        """Initialize KmerScorer object.

        """
        # self.score_map = {}  # only necessary if score kmers != vec kmers
        self.probabilities = []
        self.kmers = []

    def load(self, score_file, family):
        """Load probabilities from score file.

        Parameters
        ----------
        score_file : str
            /path/to/score_file in CSV format.
        family : str
            Family name from which to load scores.

        """
        scores = pd.read_csv(score_file)
        scores = scores[scores['label'] == family]
        # self.score_map = {
        #     kmer: score for kmer, score
        #     in zip(scores['kmer'], scores['score'])
        # }
        self.probabilities = scores['score'].values#.reshape(-1, 1)
        self.kmers = scores['kmer'].values

    def score(self, vec, vec_kmers=None, scaler=False, **scaler_kwargs):
        """Score vectors using loaded probability scores.

        Parameters
        ----------
        vec : numpy.ndarray
            Array of shape (p, q) of p vectors containing q kmers.
        vec_kmers : list of str or numpy.ndarray of str
            Array of shape (q,) containing q kmers.
        scaler : bool (default: False)
            If True, performs scaling to reduce kmer features.
        **scaler_kwargs : dict
            Keyword arguments to pass to :class: `~KmerScoreScaler()`.

        Returns
        -------
        numpy.ndarray
            Array of shape (p,) containing p scores.

        """
        # note that `vec` must be an array of arrays

        # this check doesn't work for single vec... why? raise err?
        # if np.array(vec).shape[1] != len(self.score_matrix):
        #     raise ValueError(
        #         "Kmer vector shape mismatch:"
        #         f" basis set shape is ({len(self.score_matrix)}, ) but"
        #         f" input vector shape is {np.array(vec).shape}."
        #     )

        # check that score labels match vec labels
        if vec_kmers is not None:
            if any(
                [
                    vec_kmers[i] != self.kmers[i]
                    for i in range(len(vec_kmers))
                ]
            ):
            # for i, label in enumerate(vec_labels):
                raise ValueError(
                    "Label order mismatch - `self.kmers` order is"
                    " not the same as `vec_labels`."
                )

        # to do: change to numpy matrix operations
        # note that the scaler is invariant w/ vec
        total_score = apply_feature_probabilities(
            vec, self.probabilities, scaler=scaler, **scaler_kwargs
        )
        return total_score


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


def _get_kmer_presence(feature_matrix, label_match):
    """Compute kmer presence of a given label in sequences.

    Parameters
    ----------
    feature_matrix : array (list, numpy.ndarray, pandas.Series)
        Feature matrix of shape (n, m), where each row represents
        a kmer and each column represents a sequence.
        In other words, m must equal len(labels).
    label_match : list or numpy.ndarray or pandas.Series
        Array of shape (m, ) describing label matches.
        e.g. contains 1 when kmer matches a label, and 0 otherwise

    Returns
    -------
    numpy.ndarray
        Array describing kmer presence of a given label

    """
    # convert feature count matrix into binary presence/absence
    feature_matrix = _binarize(feature_matrix)

    # for each row in feature matrix, compute each value's probability
    # no multiprocessing-- using numpy native parallelization
    result = list()
    for kmer in feature_matrix:
        kmer_match = np.equal(kmer, 1)  # array of presence matches
        ones = np.ones(len(kmer))
        match = np.multiply(np.multiply(kmer_match, label_match), ones)
        result.append(np.sum(match))
    return result


def _score(p_pos, p_neg, w=1.0):
    """Score probability vs. weighted negative probability.

    e.g. P(pos) - (w * P(neg))

    Parameters
    ----------
    p_pos : float
        Probability of positive assignment.
    p_neg : float
        Probability of negative assignment.
    w : float (default: 1.0)
        Weight for negative assignment probability.

    Returns
    -------
    float
        Total probability of positive identification, minus
        contribution to negative identification.

    """
    # p_in = results[in_label][col]
    # p_out = np.sum([results[fam]['probability'] for fam in out_labels], axis=0)
    # p_bg = bg_labels
    return p_pos - (w * p_neg)  # - (w_bg * p_bg)


def score(results, in_label, out_labels, w=1.0, col='probability', neg_only=False):
    """Score kmer from kmer family probabilities.

    The scoring method used is as follows:
        Score = P(in family) - (w_out * P(out of family)) - P(bg)

    Parameters
    ----------
    results : dict or pandas.DataFrame
        Prior results containing all kmer family probabilities.
        Requires nested dict or multi-index, e.g.:
            results[family][`col`] = pos kmer family probability
    in_label : str
        Name of "in-family" label.
    out_labels : list or array-like
        Array of "out-of-family" labels.
    w : float (default: 1.0)
        Multiplier for out-of-family probability subtraction.
    col : str (default: 'probability')
        Name of column containing kmer probabilities.
    neg_only : bool (default: False)
        If True, only returns the negative probability (-w * P_out).

    Returns
    -------
    float
        Score for in-family kmer assignment.

    """
    p_in = results[in_label][col]
    p_out = np.sum([results[fam]['probability'] for fam in out_labels], axis=0)

    # negative probability only
    if neg_only:
        return _score(0, p_out, w=w)

    return _score(p_in, p_out, w=w)


def _parse_score_method(method, bg=None, **kwargs):
    # restrict method options
    if method not in ["default", "bg_only", "both"]:
        raise ValueError("Parameter `method` must be one of the"
                         " following: 'default', 'bg_only', 'both'.")

    # check that background is provided when required
    if (method in ["bg_only", "both"]) and (bg is None):
        raise ValueError("Background sequence vectors must be"
                         " provided if `method`='bg_only' or"
                         " 'both'.")

    # define score methods
    method2name = {"default": "default_score",
                   "bg": "bg_subtracted_score",
                   "both": ["default_score", "bg_subtracted_score", "combined_score"]}

    # parse scoring methods
    if method == "default":
        # compute score based on (label, not label) assignment
        # weight = 1 / (len(labels) - 1)
        # for l in labels:
        #     o = labels[labels != l]  # other labels
        #     results[l]['score'] = score(results, l, o, w=weight)
        return "default_score", score
    elif method == "bg_only":
        return "bg_subtracted_score", score
    return "combined_score", score


def feature_class_probabilities(feature_matrix, labels,
                                kmers=None, bg_matrix=None,
                                method="default"):
    """Calculate probabilities for features being in a defined class.

    Note: only coded to work for the binary case (2 classes).

    Parameters
    ----------
    feature_matrix : array (list, numpy.ndarray, pandas.Series)
        Feature matrix of shape (n, m), where each row represents
        a kmer and each column represents a sequence.
        In other words, m must equal len(labels).
    labels : list or numpy.ndarray or pandas.Series
        Class labels of shape (m, ) describing feature matrix.
        Must have as many entries as the number of feature columns.
    kmers : list or None (default: None)
        Optional kmer identifiers mapping to kmer columns (shape (n,)).
    method : str (default: "default")
        Specify scoring method:
            "default" : Score in-family vs. out-of-family assignment
                P(family) - (weight * P(out of family))
            "bg_only" : Score in-family vs. background
                P(family) - P(background)
            "both" : Score in-family vs. out-of-family and background
                P(family) - P(background) - (weight * P(out of family))
    bg_matrix : array or None (default: None)
        Feature matrix of shape (n, o), where each row represents
        a kmer and each column represents a background sequence.

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

    # parse score function
    method_name, score_function = _parse_score_method(method)

    # kmer labels
    if kmers is None:
        kmers = np.array([n for n in range(len(feature_matrix))])
    if len(kmers) != len(feature_matrix):
        # print(len(kmers), len(feature_matrix), feature_matrix.shape)
        raise ValueError("Kmer array shape is mismatched.")

    # check that labels are the same size as number of examples
    if len(feature_matrix.T) != len(labels):
        raise ValueError("Input shapes are mismatched.")

    # get only unique labels
    unique_labels = np.unique(labels)

    # iteratively score all feature in input matrix
    n_seqs = [len(np.hstack(np.where(labels == lname))) for lname in unique_labels]
    results = {label_name: {} for label_name in labels}
    for i, l in enumerate(unique_labels):
        label_match = np.equal(labels, l)

        # get score arrays based on label match, normalized by #seqs
        presence = _get_kmer_presence(feature_matrix, label_match)
        norms = (1 / n_seqs[i]) * np.ones(len(presence))

        # if no kmers specified, save numerical range
        results[l]['kmer'] = np.array(kmers)
        results[l]['count'] = np.asarray(presence, dtype=int)
        results[l]['probability'] = np.asarray(presence * norms, dtype=float)

    # compute score based on (label, not label) assignment
    label_weight = len(unique_labels)  # if 1 label, no weight (w=1)
    if len(unique_labels) > 1:  # if >1 label, weight = 1 / (n - 1)
        label_weight = label_weight - 1
    weight = 1 / label_weight
    for l in unique_labels:
        o = unique_labels[unique_labels != l]  # other labels
        results[l]['score'] = score(results, l, o, w=weight)

    # reformat as long-form dataframe
    results = pd.DataFrame(results).T.reset_index().rename(
        columns={'index': 'label'}
        )
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
        Keyword arguments for KmerScoreScaler().

    Returns
    -------
    numpy.ndarray
        Feature matrix weighted

    """
    # optionally apply scaling
    if scaler:
        scaler = KmerScoreScaler(**kwargs)
        scaler.fit(scores)
        scores = scaler.transform(scores)
        feature_matrix = scaler.transform(feature_matrix)

    feature_matrix = _binarize(feature_matrix)
    score_matrix = scores * feature_matrix
    score_totals = np.array([np.sum(arr) for arr in score_matrix])

    return score_totals


def append_cv_results(
    data, scores, method, family=None, alphabet=None, config=None
):
    # collate ROC-AUC results
    data['family'] += [family] * config['model']['cv']
    data['alphabet_name'] += [alphabet.lower()] * config['model']['cv']
    data['k'] += [config['k']] * config['model']['cv']
    data['scoring'] += [method] * config['model']['cv']
    data['score'] += score
    # data['score'] += list(cross_val_score(clf, X, y, cv=cv, scoring='roc_auc'))
    data['cv_split'] += [i + 1 for i in range(config['model']['cv'])]
    return data
