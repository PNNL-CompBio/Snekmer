"""score: Kmer-based scoring and similarity analysis.

author: @christinehc, @biodataganache

"""
# imports
import numpy as np
import pandas as pd

from .vectorize import KmerBasis
from .utils import to_feature_matrix
from sklearn.metrics.pairwise import pairwise_distances


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

    def __init__(self, n=100):
        """Initialize KmerScaler object.

        """
        self.n = n
        self.input_shape = None
        self.scores = None
        self.basis_index = None
        self.basis_score = dict()

    def fit(self, scores):
        """Fit scaler based on list of scores.

        Parameters
        ----------
        scores : list or numpy.ndarray
            Array of scores for kmer features.
        threshold : float
            Numerical cutoff for kmer feature scores.

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
            raise ValueError("One of either `threshold` or `n` must" " be specified.")

        # store basis set and indices as attributes
        self.basis_index = indices
        self.basis_score = [
            (i, score) for i, score in zip(self.basis_index, scores[self.basis_index])
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
        if self.basis_index is None or not any(self.basis_index):
            raise AttributeError(
                "Kmer basis set not defined;" " must fit scores to scaler."
            )
        # if not isinstance(array, np.ndarray):
        #     array = np.array(array)
        if isinstance(array[0], list):
            array = to_feature_matrix(array)

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


def score(results, in_label, out_labels, w=1.0, col="probability", neg_only=False):
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
    p_out = np.sum([results[fam]["probability"] for fam in out_labels], axis=0)

    # negative probability only
    if neg_only:
        return _score(0, p_out, w=w)

    return _score(p_in, p_out, w=w)


def _parse_score_method(method, bg=None, **kwargs):
    """Translate method name into corresponding score method(s).

    Parameters
    ----------
    method : str
        Method name:
            'default' - Default score method
            'bg_only' - Score sequences as though they are all background
            'both' - Compute both score types
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
    if method not in ["default", "bg_only", "both"]:
        raise ValueError(
            "Parameter `method` must be one of the"
            " following: 'default', 'bg_only', 'both'."
        )

    # check that background is provided when required
    if (method in ["bg_only", "both"]) and (bg is None):
        raise ValueError(
            "Background sequence vectors must be"
            " provided if `method`='bg_only' or"
            " 'both'."
        )

    # define score methods
    method2name = {
        "default": "default_score",
        "bg_only": "bg_subtracted_score",
        "both": ["default_score", "bg_subtracted_score", "combined_score"],
    }

    # parse scoring methods
    if method == "both":
        return "combined_score", score
        # compute score based on (label, not label) assignment
        # weight = 1 / (len(labels) - 1)
        # for l in labels:
        #     o = labels[labels != l]  # other labels
        #     results[l]['score'] = score(results, l, o, w=weight)
    return method2name[method], score


def feature_class_probabilities(
    feature_matrix, labels, kmers=None, bg_matrix=None, method="default"
):
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
        results[l]["kmer"] = np.array(kmers)
        results[l]["count"] = np.asarray(presence, dtype=int)
        results[l]["probability"] = np.asarray(presence * norms, dtype=float)

    # compute score based on (label, not label) assignment
    label_weight = len(unique_labels)  # if 1 label, no weight (w=1)
    if len(unique_labels) > 1:  # if >1 label, weight = 1 / (n - 1)
        label_weight = label_weight - 1
    weight = 1 / label_weight
    for l in unique_labels:
        o = unique_labels[unique_labels != l]  # other labels
        results[l]["score"] = score(results, l, o, w=weight)

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
    score_totals = np.array([np.sum(arr) for arr in score_matrix])

    return score_totals


class KmerScorer:
    """Score kmer vectors based on summed probability scores.

    Attributes
    ----------
    kmers : :obj:`~KmerBasis`
        Basis set object for kmer standardization.
    matrix : numpy.ndarray
        Feature matrix for all sequences.
    labels : list or array-like of str
        Array of all labels in overall dataset.
    primary_label : str
        Primary label identifier for fitted sequence data.
    i_background : dict of list
        Description of attribute `i_background`.
    i_label : dict
        Description of attribute `i_label`.
    scaler : type
        Description of attribute `scaler`.
    scores : type
        Description of attribute `scores`.
    score_norm : type
        Description of attribute `score_norm`.
    probabilities : type
        Description of attribute `probabilities`.

    """

    def __init__(self):
        self.kmers = KmerBasis()
        self.primary_label = None
        self.scaler = None
        self.scores = {"sample": {}, "background": {}}
        self.score_norm = None
        self.probabilities = {"sample": {}, "background": {}}

    # load list of kmers and both seq and bg feature matrices
    def fit(
        self,
        kmers,
        data,
        label,
        label_col="family",
        bg_col="background",
        vec_col="vector",
        **scaler_kwargs,
    ):
        """Short summary.

        Parameters
        ----------
        kmers : type
            Description of parameter `kmers`.
        data : type
            Description of parameter `data`.
        label : type
            Description of parameter `label`.
        label_col : type
            Description of parameter `label_col`.
        bg_col : type
            Description of parameter `bg_col`.
        **scaler_kwargs : dict
            Keyword arguments for
            :obj:`~snekmer.model.KmerScoreScaler`.

        Returns
        -------
        type
            Description of returned object.

        """
        # save primary family label
        self.primary_label = label

        # step 00: set kmer basis
        self.kmers.set_basis(kmers)

        # step 0: get indices of sample and background sequences
        i_background = {
            "sample": list(data.index[~data[bg_col]]),
            "background": list(data.index[data[bg_col]]),
        }


        # step 0: get feature matrix and all labels
        labels = data[label_col].values
        matrix = to_feature_matrix(data[vec_col])
        #print(data[vec_col][1])

        #matrix,labels = to_feature_matrix(data[vec_col])
        print(matrix.shape)

        # step 0: get indices for label (family) ids
        # i_label = {
        #     ll: list(data.index[data[label_col] == ll])
        #     for ll in data[label_col].unique()
        # }

        # step 1: score sample sequences and fit score scaler
        sample_matrix = matrix[i_background["sample"]]
        s_labels = labels[i_background["sample"]]

        probas = feature_class_probabilities(
            sample_matrix.T, s_labels, kmers=self.kmers.basis
        )
        self.probabilities["sample"] = probas[probas["label"] == self.primary_label][
            "score"
        ].values

        # step 2: fit scaler to the sample data (ignore the background)
        self.scaler = KmerScoreScaler(**scaler_kwargs)
        self.scaler.fit(
            probas[probas["label"] == label]["probability"]
        )  # probas is an ordered list of scores from kmers. returns indices for these scores

        # step 3: compute background
        background_matrix = matrix[i_background["background"]]
        bg_labels = labels[i_background["background"]]
        if len(background_matrix) > 0:
            bg_probas = feature_class_probabilities(
                background_matrix.T, bg_labels, kmers=self.kmers.basis
            )
            self.probabilities["background"] = bg_probas[
                bg_probas["label"] == self.primary_label
            ]["score"].values

            # step 3.1: background family probability scores
            for bl in np.unique(bg_labels):
                bg_scores = bg_probas[bg_probas["label"] == bl]["score"].values

                # normalize by max bg score
                # bg_norm = np.max(bg_scores)

                # get background scores for sequences
                bg_only_scores = apply_feature_probabilities(
                    matrix, bg_scores, scaler=self.scaler
                )

                # normalize by max bg score attained by a sequence
                bg_norm = np.max(bg_only_scores)

                self.scores["background"][bl] = bg_only_scores / bg_norm

        # step 3.2: assign family probability scores
        for sl in np.unique(s_labels):
            scores = probas[probas["label"] == sl]["score"].values

            # normalize by sum of all positive scores
            # norm = np.sum([s for s in scores if s > 0])

            # include background sequences for score generation
            total_scores = apply_feature_probabilities(
                matrix, scores, scaler=self.scaler
            )

            # normalize by max score
            norm = np.max(total_scores)
            norm = max(norm, 1.0)  # prevent score explosion?
            if sl == self.primary_label:
                self.score_norm = norm

            # assign percent score based on max positive score
            self.scores["sample"][f"{sl}_score"] = total_scores / norm

            # weight family score by (1 - normalized bg score)
            if sl in np.unique(bg_labels):
                self.scores["sample"][f"{sl}_score_background_weighted"] = [
                    total * (1 - bg)
                    for total, bg in zip(
                        self.scores["sample"][f"{sl}_score"],
                        self.scores["background"][sl],
                    )
                ]

                # old scoring method
                self.scores["sample"][f"{sl}_background_subtracted_score"] = [
                    total - bg
                    for total, bg in zip(
                        self.scores["sample"][f"{sl}_score"],
                        self.scores["background"][sl],
                    )
                ]

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
        array = self.kmers.transform(array, kmers)
        return (
            apply_feature_probabilities(
                array, self.probabilities["sample"], scaler=self.scaler
            )
            / self.score_norm
        )
