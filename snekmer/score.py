"""score: Kmer-based scoring and similarity analysis.

author: @christinehc, @biodataganache

"""
# imports
from typing import Optional, Union
import numpy as np
import pandas as pd

from ._version import __version__
from .vectorize import KmerBasis
from numpy.typing import ArrayLike, NDArray
from sklearn.metrics.pairwise import pairwise_distances

# define variables
SCORE_METHODS = ["background", "bg", "family", "f", "combined"]
METHOD2NAME = {
    "background": "bg_subtract",
    "bg": "bg_subtract",
    "family": "family_subtract",
    "f": "family_subtract",
    "combined": "bg_family_subtract",
}


# scaling object for reducing kmers used in scoring
class KmerScoreScaler:
    """Scale and reduce kmer set based on scores.

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

    def __init__(self, n: int = 100):
        """Initialize KmerScaler object."""
        self.n = n
        self.input_shape = None
        self.scores = None
        self.basis_index = None
        self.basis_score = dict()

    def fit(self, scores: ArrayLike):
        """Fit scaler based on list of scores.

        Parameters
        ----------
        scores : ArrayLike
            Array of scores for kmer features.

        Returns
        -------
        None
            Fits KmerScoreScaler() object.

        """
        if not isinstance(scores, np.ndarray):
            scores = np.array(scores)
        self.scores = scores
        self.input_shape = scores.shape

        # option 1: set score threshold as limit -- REMOVED
        # if threshold is not None:
        #     indices = np.where(np.array(scores > threshold))[0]

        if self.n is not None:
            # option 2: set # of features as limit
            if (isinstance(self.n, int)) and (self.n >= 1):
                indices = scores.ravel().argsort()[: -self.n - 1 : -1]
            # option 3: set % of features as limit
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
            raise ValueError("Feature cutoff `n` must" " be specified.")

        # store basis set and indices as attributes
        self.basis_index = indices
        self.basis_score = [
            (i, score) for i, score in zip(self.basis_index, scores[self.basis_index])
        ]
        return

    def transform(self, array: ArrayLike) -> ArrayLike:
        """Reduce vector to kmer basis set.

        Parameters
        ----------
        array : ArrayLike
            Unscaled kmer vector or matrix.

        Returns
        -------
        ArrayLike
            Scaled kmer vector or matrix.

        """
        if self.basis_index is None or not any(self.basis_index):
            raise AttributeError(
                "Kmer basis set not defined;" " must fit scores to scaler."
            )

        # enforce array shape
        array = np.asarray(array)
        array = np.vstack(array)

        # input vec must match original input shape for indexing
        if len(array.shape) < 2:  # account for 1D arrays
            array_size = len(array)
            input_size = self.input_shape[0]
            if array_size != input_size:
                raise ValueError(
                    f"Input vector shape {array.shape}"
                    f" does not match fitted vector shape"
                    f" {self.input_shape}."
                )
            return array[self.basis_index]

        # ND arrays
        array_size = array.shape
        input_size = self.input_shape
        return array[:, self.basis_index]


# functions
def connection_matrix_from_features(
    feature_matrix: ArrayLike, metric: str = "jaccard"
) -> NDArray:
    """Calculate similarities based on features output from main().

    Parameters
    ----------
    feature_matrix : ArrayLike
        Feature matrix in the form of rows (sequences), columns
        (kmers), where values are kmer counts for each protein.
    metric : str
        Metric used for pairwise distance computation (see
            sklearn.metrics.pairwise documentation).

    Returns
    -------
    NDArray
        Square matrix with the similarity scores of pairwise
        relationships between proteins.

    """
    if metric == "jaccard":
        sim = 1 - pairwise_distances(feature_matrix, metric="hamming")  # .T?
    else:
        sim = pairwise_distances(feature_matrix, metric=metric)
    return sim


def _apply_probability(kmer: int, label: str, compare_label: str) -> int:
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


def _binarize(feature_matrix: ArrayLike) -> ArrayLike:
    """Convert feature counts into binary presence/absence.

    Parameters
    ----------
    feature_matrix : ArrayLike
        Feature matrix, where each row represents a kmer and each
        column represents a sequence

    Returns
    -------
    ArrayLike
        Binarized presence/absence matrix.

    """
    return (feature_matrix > 0) * 1.0


def _get_kmer_presence(feature_matrix: ArrayLike, label_match: ArrayLike) -> NDArray:
    """Get kmer presence in full feature matrix for one family.

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

    # for each kmer, compute presence in seqs
    # then account for seqs matching a given label
    # end result is presence weighted by family membership
    result = list()
    for kmer in feature_matrix:
        kmer_match = np.equal(kmer, 1)  # boolean array (present/not)
        ones = np.ones(len(kmer))
        match = np.multiply(np.multiply(kmer_match, label_match), ones)
        result.append(np.sum(match))
    return result


def _score(p_pos: float, p_neg: float, w: float = 1.0) -> float:
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


def score_old(results, in_label, out_labels, w=1.0, col="probability"):
    """DEPRECATED: Score kmer from kmer family probabilities.

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

    Returns
    -------
    float
        Score for in-family kmer assignment.

    """
    # assume family wt is second list element unless list of len 1
    if isinstance(w, list) and len(w) == 1:
        w = w[0]
    if isinstance(w, list) and len(w) > 1:
        w = w[1]
    p_in = results[in_label][col]
    p_out = np.sum([results[fam][col] for fam in out_labels], axis=0)

    return _score(p_in, p_out, w=w)


def score(
    results: Union[dict, pd.DataFrame],
    in_label: str,
    out_labels: ArrayLike,
    background: Optional[NDArray] = None,
    w_bg: float = 0.25,
    w_out: float = 1.0,
    col: str = "probability",
    method: str = "background",
) -> float:
    """Score kmer using specified scoring method.

    Parameters
    ----------
    results : Union[dict, pd.DataFrame]
        Prior results containing all kmer family probabilities.
        Requires nested dict or multi-index, e.g.:
            results[family][`col`] = pos kmer family probability
    in_label : str
        Name of "in-family" label.
    out_labels : ArrayLike
        Array of "out-of-family" labels.
    background : Optional[NDArray], optional
        Background matrix, by default None
    w_bg : float, optional
        Weight multiplier for background subtraction,
            by default 0.25
    w_out : float, optional
        Weight multiplier for family subtraction, by default 1.0
    col : str, optional
        Name of column containing kmer probabilities,
            by default "probability"
    method : str, optional
        Scoring method, by default "background"
            Accepts ["combined", "background", "bg", "family", "f"]

    Returns
    -------
    float
        _description_

    Raises
    ------
    ValueError
        _description_
    """
    p_in = results[in_label][col]
    p_out = np.sum([results[fam][col] for fam in out_labels], axis=0)

    if method in ["combined", "bg", "background"]:
        bg = np.abs(_score(0, background, w=w_bg))
        if background is None:
            raise ValueError(
                f"Background matrix `background` must be provided for `method`='{method}'."
            )

    if method in ["combined", "family", "f"]:
        out = abs(_score(0, p_out, w=w_out))

    if method == "combined":
        return p_in - out - bg

    if method in ["family", "f"]:
        return p_in - out

    return p_in - bg


def _parse_score_method(method, bg=None, **kwargs):
    """Translate method name into corresponding score method(s).

    Parameters
    ----------
    method : str
        Method name:
            'background'/'bg' - Score with background subtracted
            'family'/'f' - Score with (in, out) family labels
                           (i.e. out-of-family (OOF) subtracted)
            'combined' - Subtract both background and OOF
    bg : type
        Description of parameter `bg`.
    **kwargs : type
        Description of parameter `**kwargs`.

    Returns
    -------
    type
        Description of returned object.

    """
    # restrict method options
    if method not in SCORE_METHODS:
        raise ValueError(
            f"Parameter `method` must be one of the following: {SCORE_METHODS}."
        )

    # check that background is provided when required
    if (method in ["background", "bg", "combined"]) and (bg is None):
        raise ValueError(
            "Background sequence vectors must be" f" provided for `method`='{method}'."
        )

    # parse scoring methods
    return METHOD2NAME[method], score


def feature_class_probabilities(
    feature_matrix: ArrayLike,
    labels: ArrayLike,
    kmers: Optional[list] = None,
    bg: Optional[ArrayLike] = None,
    method: str = "background",
    weight_bg: float = 0.25,
):
    """Calculate probabilities for features being in a defined class.

    Note: only coded to work for the binary case (2 classes).

    Parameters
    ----------
    feature_matrix : ArrayLike
        Feature matrix of shape (n, m), where each row represents
        a kmer and each column represents a sequence.
        In other words, m must equal len(labels).
    labels : ArrayLike
        Class labels of shape (m, ) describing feature matrix.
        Must have as many entries as the number of feature columns.
    kmers : list, optional
        Optional kmer identifiers mapping to kmer columns (shape (n,)),
            by default None
    bg : ArrayLike, optional (default: None)
        Feature matrix of shape (n, o), where each row represents
        a kmer and each column represents a background sequence.
    method : str (default: "default")
        Specify scoring method:
            "background"/"bg" : Score in-family vs. background
                P(family) - P(background)
            "family"/"f": Score in-family vs. out-of-family assignment
                P(family) - (weight * P(out of family))
            "combined" : Score in-family vs. out-of-family and background
                P(family) - P(background) - (weight * P(out of family))
    weight_bg : float, optional
        Weight multiplier for background subtraction,
            by default 0.25


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
    # method_name, score_function, bg_matrix = _parse_score_method(method, bg=bg)

    # kmer labels
    kmers = np.asarray(kmers)
    if kmers is None:
        kmers = np.arange(len(feature_matrix))
    if len(kmers) != len(feature_matrix):
        raise ValueError(
            "Input kmer array shape does not match kmer basis"
            f" ({len(feature_matrix)} vs. {len(kmers)})."
        )

    # check that labels are the same size as number of examples
    if len(feature_matrix.T) != len(labels):
        raise ValueError(
            "Input kmer array size does not match labels"
            f" ({len(feature_matrix)} vs. {len(labels)})."
        )

    # if bg matrix exists, verify shape fits no. kmers
    if bg is not None:
        if len(bg.T) != len(kmers):
            raise ValueError(
                "Background matrix shape does not match kmer basis"
                f" ({len(bg.T)} vs. {len(kmers)})."
            )

    # get only unique labels
    unique_labels = np.unique(labels)

    # iteratively score all feature in input matrix
    n_seqs = [np.sum(labels == lname) for lname in unique_labels]  # count seqs /label
    results = {label_name: {} for label_name in labels}
    for i, l in enumerate(unique_labels):
        label_match = np.equal(labels, l)

        # get score arrays based on label match, normalized by #seqs
        presence = _get_kmer_presence(feature_matrix, label_match)
        norms = (1 / n_seqs[i]) * np.ones(len(presence))

        # if no kmers specified, save numerical range
        if bg is not None:
            background = _get_kmer_presence(bg, np.ones(len(bg)))
            bg_norms = (1 / np.sum(background)) * np.ones(len(background))
            p_background = np.asarray(background * bg_norms, dtype=float)
        else:
            p_background = None

        results[l]["kmer"] = kmers
        results[l]["count"] = np.asarray(presence, dtype=int)
        results[l]["probability"] = np.asarray(presence * norms, dtype=float)

    # compute score based on (label, not label) assignment
    label_weight = len(unique_labels)  # if 1 label, no weight (w=1)
    if len(unique_labels) > 1:  # if >1 label, weight = 1 / (n - 1)
        label_weight = label_weight - 1
    weight_out = 1 / label_weight
    for l in unique_labels:
        o = unique_labels[unique_labels != l]  # other labels
        results[l][f"score_{METHOD2NAME[method]}"] = score(
            results,
            l,
            o,
            method=method,
            background=p_background,
            w_bg=weight_bg,
            w_out=weight_out,
        )

    # reformat as long-form dataframe
    results = pd.DataFrame(results).T.reset_index().rename(columns={"index": "label"})
    results = results.set_index(["label"]).apply(pd.Series.explode).reset_index()
    return results


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
        scaler = KmerScoreScaler(**kwargs)
        scaler.fit(scores)
    if isinstance(scaler, KmerScoreScaler):
        # scaler.fit(scores)
        scores = scaler.transform(scores)
        feature_matrix = scaler.transform(feature_matrix)

    feature_matrix = _binarize(feature_matrix)

    score_matrix = scores * feature_matrix
    score_totals = np.sum(score_matrix, axis=1)

    return score_totals


# new KmerScorer obj
class KmerScorer:
    """Score kmer vectors based on summed probability scores.

    Attributes
    ----------
    snekmer_version : str
        Snekmer version used to generate current object.
    kmers : :obj:`~KmerBasis`
        Basis set object for kmer standardization.
    label : str
        Primary label identifier for fitted sequence data.
    scaler : :obj:`~KmerScoreScaler`
        Scaler object to reduce kmer set based on scores.
    probabilities : ArrayLike
        Kmer probabilities calculated from freqs.
    score_norm : float
        Normalization constant for scores.
    scores : pandas.DataFrame
        Output scores.

    """

    def __init__(self, method="background"):
        self.snekmer_version = __version__
        self.method = method
        self.kmers = KmerBasis()
        self.label = None
        self.scaler = None
        self.probabilities = None
        self.score_norm = None
        self.scores = None
        self.score_col = f"score_{METHOD2NAME[self.method]}"

    # load list of kmers and both seq and bg feature matrices
    def fit(
        self,
        kmers: ArrayLike,
        data: pd.DataFrame,
        label: str,
        bg: Optional[ArrayLike] = None,
        weight_bg=0.25,
        label_col: str = "family",
        vec_col: str = "vector",
        **scaler_kwargs,
    ) -> None:
        """Fit and score data.

        Parameters
        ----------
        kmers : ArrayLike
            _description_
        data : pd.DataFrame
            _description_
        label : str
            _description_
        label_col : str, optional
            _description_, by default "family"
        bg : ArrayLike, optional
            _description_, by default None
        vec_col : str, optional
            _description_, by default "vector"
        """
        # save primary family label
        self.label = label

        # step 1: use kmer set to define basis
        self.kmers.set_basis(kmers)

        # step 2: get feature matrix and all labels
        labels = data[label_col].values
        x = len(data[vec_col])
        y = len(data[vec_col][0])

        matrix = np.zeros(x * y).reshape((x, y))
        for i in range(x):
            for j in range(y):
                value = data[vec_col][i]
                value = value[j]
                matrix[i, j] = value
        probas = feature_class_probabilities(
            matrix.T,
            labels,
            kmers=self.kmers.basis,
            method=self.method,
            bg=bg,
            weight_bg=weight_bg,
        )
        self.probabilities = probas[probas["label"] == self.label][
            self.score_col
        ].to_numpy()

        # step 3: fit scaler to the sample data (ignore the background)
        self.scaler = KmerScoreScaler(**scaler_kwargs)
        self.scaler.fit(
            probas[probas["label"] == label]["probability"]
        )  # probas is an ordered list of scores from kmers. returns indices for these scores

        # get probability scores for each label
        scores = apply_feature_probabilities(
            matrix, self.probabilities, scaler=self.scaler
        )

        # normalize by sum of all positive scores
        # norm = np.sum([s for s in scores if s > 0])

        # normalize by max score
        self.score_norm = max(np.max(scores), 1.0)  # prevent score explosion?

        # assign percent score based on max positive score
        self.scores = np.array(scores) / self.score_norm

    # score new input sequences
    def predict(self, array, kmers):
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
        array = self.kmers.transform(array, kmers)
        return (
            apply_feature_probabilities(
                array,
                self.probabilities,
                scaler=False,  # self.scaler
            )
            / self.score_norm
        )
