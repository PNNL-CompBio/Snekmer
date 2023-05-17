"""score: Kmer-based scoring and similarity analysis.

author: @christinehc, @biodataganache

"""
# imports
import numpy as np
import pandas as pd

from ._version import __version__

# from .vectorize import KmerBasis
from .utils import to_feature_matrix
from sklearn.metrics.pairwise import pairwise_distances


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
    return np.array(result)


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


def score(results, in_label, labels, w=1.0, col="probability", neg_only=False):
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
    labels : list or array-like
        Array all labels.
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
    out_labels = np.unique([l for l in labels if l != in_label])
    p_in = results[in_label][col]
    p_out = np.sum([results[fam]["probability"] for fam in out_labels], axis=0)

    # negative probability only
    if neg_only:
        return _score(0, p_out, w=w)

    return _score(p_in, p_out, w=w)


# from Jason
# This code will apply the feature class probabilities derived in the previous functions
#   to a new set of examples as a weighted score. The idea here is to use a very simple
#   approach to classification.
def apply_feature_probabilities(feature_matrix, scores, scaler=False, **kwargs):
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
    scaler : obj or bool (default: False)
        If True, performs scaling to reduce kmer features.
        If False, does not perform scaling to reduce kmer features.
        If obj is supplied, the scaling obj is applied to reduce
            kmer features.
    **kwargs : dict
        Keyword arguments for KmerScoreScaler().
    Returns
    -------
    numpy.ndarray
        Feature matrix weighted
    """
    # optionally apply scaling
    if scaler and isinstance(scaler, bool):
        scaler = KmerSetReducer(**kwargs)
        scaler.fit(scores)
    if isinstance(scaler, KmerSetReducer):
        # scaler.fit(scores)
        scores = scaler.transform(scores)
        feature_matrix = scaler.transform(feature_matrix)

    feature_matrix = _binarize(feature_matrix)

    score_matrix = scores * feature_matrix
    score_totals = np.array([np.sum(arr) for arr in score_matrix])

    return score_totals


# scaling object for reducing kmers used in scoring
class KmerSetReducer:
    """Reduce kmer set based on scores/prevalence.

    Parameters
    ----------
    scores : list or numpy.ndarray
        Description of parameter `model`.
    n : int
        Feature cutoff (top n kmer features to include)

    Attributes
    ----------
    n : int or float
        Feature cutoff (top n kmer features).
        Accepts either an integer of the top n features, or a float
        indicating a percentage (n = 0.1 for the top 10% of features)
    input_shape : tuple
        Input array shape.
    kmer_basis_idx : numpy.ndarray
        List of indices from kmer feature vector making up the
        specified basis set, ordered from most to least important.
    kmer_basis_score : list of tuples
        List of tuples for kmer feature vector making up the
        specified basis set, ordered from most to least important.
        Tuples are: (index, score) for each listed feature.
    scores : numpy.ndarray
        Array of scores for kmer features.

    """

    def __init__(self, n=100):
        """Initialize KmerSetReducer object."""
        self.n = n
        self.scores = None
        self.index = None
        self.basis_score = dict()

    def fit(self, scores):
        """Fit object based on list of scores.

        Parameters
        ----------
        scores : list or numpy.ndarray
            Array of scores for kmer features.
        threshold : float
            Numerical cutoff for kmer feature scores.

        Returns
        -------
        None
            Fits KmerSetReducer() object.

        """
        if not isinstance(scores, np.ndarray):
            scores = np.array(scores)
        self.scores = scores

        # set threshold at which to keep kmers
        if self.n is not None:
            # option 1: set # of features as limit
            if (isinstance(self.n, int)) and (self.n >= 1):
                indices = scores.ravel().argsort()[: -self.n - 1 : -1]
            # option 2: set % of features as limit
            elif (isinstance(self.n, float)) and (self.n < 1):
                n = int(np.floor(self.n * len(scores)))
                indices = scores.ravel().argsort()[: -n - 1 : -1]
            else:
                raise ValueError(
                    "Invalid input format for `n`"
                    " (must be either int > 1, or"
                    " 0.0 < float < 1.0)."
                )
        else:
            raise ValueError("One of either `threshold` or `n` must be specified.")

        # store basis set and indices as attributes
        self.index = sorted(indices)
        self.basis_score = [
            (i, score) for i, score in zip(self.index, scores[self.index])
        ]
        return

    def transform(self, array):
        """Reduce vector to kmer basis set.

        Parameters
        ----------
        array : list or numpy.ndarray
            Unscaled kmer vector or matrix.

        Returns
        -------
        numpy.ndarray
            Scaled kmer vector or matrix.

        """
        if self.index is None or not any(self.index):
            raise AttributeError(
                "Kmer basis set not defined; must fit scores to scaler."
            )
        # if not isinstance(array, np.ndarray):
        #     array = np.array(array)
        if isinstance(array[0], list):
            array = to_feature_matrix(array)

        # input vec must match original input shape for indexing
        if len(array.shape) < 2:  # account for 1D arrays
            array_size = len(array)
            input_size = self.scores.shape[0]
            if array_size != input_size:
                raise ValueError(
                    f"Input vector shape {array.shape}"
                    f" does not match fitted vector shape"
                    f" {self.scores.shape}."
                )
            return array[self.index]

        # ND arrays
        array_size = array.shape
        input_size = self.scores.shape
        return array[:, self.index]


# assume matrix is already pre-harmonized
class KmerScorer:
    """Score kmer vectors based on summed probability scores."""

    def __init__(self):
        self.label = None
        self.score_norm = None
        self.probabilities = {"sample": {}, "background": {}}
        self.snekmer_version = __version__

    # fit scorer to matrix
    def fit(self, X, y, label, codes=None):
        """_summary_

        Parameters
        ----------
        X : _type_
            _description_
        y : _type_
            _description_
        label : _type_
            _description_
        codes : _type_, optional
            _description_, by default None
        """        
        # encode labels
        # le = LabelEncoder()
        # le.fit(y)
        # Y = le.transform(y)

        # set in-family label if none specified
        # if label not in y:
        #     raise ValueError("Label not found in label array `y`.")
        self.label = label

        # TODO: account for background sequences

        # get unique labels and seq counts per label
        y_unique = np.unique(y)
        n_seqs = {
            l: len(np.hstack(np.where(y == l))) for l in y_unique
        }  # count seqs per label

        # compute score based on (label, not label) assignment
        weight = len(y_unique)  # if 1 label, no weight (w=1)
        if len(y_unique) > 1:  # if >1 label, weight = 1 / (n - 1)
            weight = weight - 1
        weight = 1 / weight

        # iteratively score all features in input matrix
        results = {l: {} for l in y_unique}
        for l in y_unique:
            # find rows that match desired label
            match = np.equal(y, l)

            # get score arrays based on label match, normalized by #seqs
            presence = _get_kmer_presence(X.T, match)
            norms = (1 / n_seqs[l]) * np.ones(len(presence))

            # if no kmers specified, save numerical range
            results[l]["kmercode"] = codes
            results[l]["count"] = np.asarray(presence, dtype=int)
            results[l]["probability"] = np.asarray(presence * norms, dtype=float)

        for l in y_unique:
            results[l]["score"] = score(results, l, y, w=weight)

            # normalize by max score
            total_score = np.array([np.sum(arr) for arr in results[l]["score"]])
            norm = np.max(total_score)
            norm = max(norm, 1.0)  # prevent score explosion?
            if l == 1:
                self.score_norm = norm

            # assign percent score based on max positive score
            results[l]["family_score"] = total_score / norm

            # weight family score by (1 - normalized bg score)
        #             if sl in np.unique(bg_labels):
        #                 self.scores["sample"][f"{sl}_score_background_weighted"] = [
        #                     total * (1 - bg)
        #                     for total, bg in zip(
        #                         self.scores["sample"][f"{sl}_score"],
        #                         self.scores["background"][sl],
        #                     )
        #                 ]
        self.score = (
            pd.DataFrame(results.values(), index=results.keys())
            .reset_index()
            .rename(columns={"index": "label"})
            .explode(["kmercode", "count", "probability", "score", "family_score"])
        )

    # score new input sequences
    def predict(self, array):
        """Predict scores of new array with separate kmer basis.

        Parameters
        ----------
        array : list or array-like
            Description of parameter `array`.
        kmers : list or array-like of str
            Ordered array of kmer strings.

        Returns
        -------
        type
            Description of returned object.

        """
        # fit new input data to the correct kmer basis
        # self.kmers.set_basis(kmers)
        # array = self.kmers.transform(array, kmers)

        # return (
        #    apply_feature_probabilities(
        #        array, self.probabilities["sample"], scaler=self.scaler
        #    )
        #    / self.score_norm
        # )
        return (
            apply_feature_probabilities(
                array,
                scorer.score[scorer.score["label"] == 1]["probability"].values,
                scaler=False,
            )
            / self.score_norm
        )
