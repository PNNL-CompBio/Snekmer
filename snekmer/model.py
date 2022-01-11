"""model: Machine learning models for Kmer family prediction.

author: @christinehc
"""
# imports
import collections.abc
import json

import numpy as np
import pandas as pd

from .io import read_output_kmers
from .utils import check_list, get_family, to_feature_matrix
from .plot import show_cv_roc_curve, show_cv_pr_curve
from .score import (
    feature_class_probabilities, apply_feature_probabilities, KmerScoreScaler
)
from .transform import KmerBasis
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier, AdaBoostClassifier
from sklearn.linear_model import LogisticRegressionCV
from sklearn.model_selection import GridSearchCV, cross_validate
from sklearn.pipeline import make_pipeline, Pipeline

# define default gridsearch param dict
MODEL_PARAMS = {
#     'scaler__n': [None, 100],
    'clf__class_weight': ['balanced'],
    'clf__min_samples_split': [2, 5, 10, 15, 20, 40],
    'clf__min_samples_leaf': [2, 5, 10, 15, 20],
    'clf__max_depth': [3, 10, 25],
  }

MODEL_NAME =  {
    'decisiontree': DecisionTreeClassifier(),
    'randomforest': RandomForestClassifier(),
    'adaboost': AdaBoostClassifier()
}


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
        self.scores = {'sample': {}, 'background': {}}
        self.score_norm = None
        self.probabilities = {'sample': {}, 'background': {}}

    # load list of kmers and both seq and bg feature matrices
    def fit(self, kmers, data, label,
            label_col='family', bg_col='background',
            **scaler_kwargs):
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
            'sample': list(data.index[~data[bg_col]]),
            'background': list(data.index[data[bg_col]])
        }

        # step 0: get feature matrix and all labels
        labels = data[label_col].values
        matrix = to_feature_matrix(data['vector'])
        # print(matrix.shape)

        # step 0: get indices for label (family) ids
        # i_label = {
        #     ll: list(data.index[data[label_col] == ll])
        #     for ll in data[label_col].unique()
        # }

        # step 1: score sample sequences and fit score scaler
        sample_matrix = matrix[i_background['sample']]
        s_labels = labels[i_background['sample']]
        probas = feature_class_probabilities(
            sample_matrix.T, s_labels, kmers=self.kmers.basis
        )
        self.probabilities['sample'] = probas[
            probas['label'] == self.primary_label
        ]['score'].values

        # step 2: fit scaler to the sample data (ignore the background)
        self.scaler = KmerScoreScaler(**scaler_kwargs)
        self.scaler.fit(
            probas[probas['label'] == label]['probability']
        )  # probas is an ordered list of scores from kmers. returns indices for these scores

        # step 3: compute background
        background_matrix = matrix[i_background['background']]
        bg_labels = labels[i_background['background']]
        if len(background_matrix) > 0:
            bg_probas = feature_class_probabilities(
                background_matrix.T, bg_labels, kmers=self.kmers.basis
            )
            self.probabilities['background'] = bg_probas[
                bg_probas['label'] == self.primary_label
            ]['score'].values

            # step 3.1: background family probability scores
            for bl in np.unique(bg_labels):
                bg_scores = bg_probas[
                    bg_probas['label'] == bl
                ]['score'].values

                # normalize by max bg score
                # bg_norm = np.max(bg_scores)

                # get background scores for sequences
                bg_only_scores = apply_feature_probabilities(
                    matrix, bg_scores, scaler=self.scaler
                )

                # normalize by max bg score attained by a sequence
                bg_norm = np.max(bg_only_scores)

                self.scores['background'][bl] = bg_only_scores / bg_norm

        # step 3.2: assign family probability scores
        for sl in np.unique(s_labels):
            scores = probas[probas['label'] == sl]['score'].values

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
            self.scores['sample'][f"{sl}_score"] = total_scores / norm

            # weight family score by (1 - normalized bg score)
            if sl in np.unique(bg_labels):
                self.scores['sample'][
                    f"{sl}_score_background_weighted"
                ] = [
                    total * (1 - bg) for total, bg in zip(
                        self.scores['sample'][f"{sl}_score"],
                        self.scores['background'][sl]
                    )
                ]

                # old scoring method
                self.scores['sample'][
                    f"{sl}_background_subtracted_score"
                ] = [
                    total - bg for total, bg in zip(
                        self.scores['sample'][f"{sl}_score"],
                        self.scores['background'][sl]
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
        return apply_feature_probabilities(
            array, self.probabilities['sample'], scaler=self.scaler
        ) / self.score_norm


# classification models for protein families
class KmerModel:
    """Classify a protein family using kmer vectors as input.

    Attributes
    ----------
    scaler :  KmerScaler object or None
        Scaler object; if None, feature vectors are not scaled
    params : dict
        Parameter dictionary for hyperparameter gridsearch.
    model : str (default: 'decisiontree')
        String identifier for classifier model.
        Mappings to Classifier object are:
            {
            'decisiontree': DecisionTreeClassifier,
            'randomforest': RandomForestClassifier,
            'adaboost': AdaBoostClassifier
            }
    step_name : str (default: 'clf')
        Optional custom name for classifier pipeline step.
    pipeline : sklearn.pipeline.Pipeline object
        Pipeline object for preprocessing and modeling.
    search : sklearn.model_selection.GridSearchCV object
        GridSearchCV object for hyperparameter searching.

    """
    def __init__(self,
                 scaler=None,
                 model='decisiontree',
                 params=MODEL_PARAMS,
                 step_name="clf"):
        """Initialize KmerModel object.

        Parameters
        ----------
        scaler :  KmerScaler object or None (default: None)
            Scaler object; if None, does not scale feature vectors
        model : Classifier object (default: DecisionTreeClassifier())
            Type of classification model.
        params : dict (default: MODEL_PARAMS)
            Parameter dictionary for hyperparameter gridsearch.
        step_name : str (default: 'clf')
            Optional custom name for classifier pipeline step.

        """
        self.scaler = scaler
        self.params = params
        self.model = MODEL_NAME[model]
        self.step_name = step_name

        self.pipeline = None
        self.search = None

    def fit(self, X, y, verbose=True):
        """Train model using gridsearch-tuned hyperparameters.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            Training input samples.
        y : array-like of shape (n_samples,) or (n_samples, n_outputs)
            Target values (class labels) as integers or strings.
        scores : array-like of shape (n_features,) or (n_features, n_outputs)
            NOT IMPLEMENTED YET-- for KmerScaler integration.
        verbose : bool (default: True)
            If True, prints best gridsearch CV results to console.

        Returns
        -------
        None
            Fits estimator.

        """
        # define pipeline and gridsearch
        self.pipeline = Pipeline(steps=[("scaler", self.scaler), (self.step_name, self.model)])
        self.search = GridSearchCV(self.pipeline, self.params)

        # use gridsearch on training data to find best parameter set
        self.search.fit(X, y)
        if verbose:
            print("best mean cross-validation score: {:.3f}".format(
                self.search.best_score_
                ))
            print("best parameters:", self.search.best_params_)

        # fit model with best gridsearch parameters
#         self.scaler.fit(scores)
        self.model.fit(X, y)
        return

    def score(self, X, y):
        """Score model on test set.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            Test samples.
        y : array-like of shape (n_samples,) or (n_samples, n_outputs)
            Predicted classes, or predict values.

        Returns
        -------
        float
            Model prediction accuracy on test set.

        """
        return self.search.score(X, y)

    def predict(self, X):
        """Assign classification based on new input vector.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            New samples.

        Returns
        -------
        int
            (0,1) classification result

        """
        return self.model.predict(X)


# functions
def format_data_df(filenames, label_name='family'):
    """Format Kmer sequence data into long-form dataframe.

    Each row of the dataframe
    Dataframe consists of the following columns:
        'seq_id': Sequence IDs
        'vec': Kmer-ized vectors
        'family' (or other label name): Pre-defined label
        '{label_0, ..., label_n}': Positive identification as a
                                   given label; 0 for label, 1
                                   for non-label. Note that this
                                   gives n separate columns.

    Parameters
    ----------
    filenames : list
        List of filenames for aggregation.
        Files must all correspond to standardized k-mer outputs.
    label_name : str
        Name of sequence label (default: "family")

    Returns
    -------
    pandas.DataFrame
        Dataframe containing formatted sequence data.

    """
    data = {'seq_id': [], 'vec': [], label_name: []}
    for fn in filenames:
        with open(fn, 'r') as f:
            tmp = json.load(f)
            data['seq_id'] += [tmp['seq_id']]
            data['vec'] += [np.array(tmp['vector'])]
            data[label_name] += [get_family(fn)]

    data[label_name] = np.hstack(
        np.array(
            [[label] * len(data['vec'][i]) for i, label in enumerate(data[label_name])],
            dtype='object'
            )
        )
    data['seq_id'] = np.hstack([np.array(s) for s in data['seq_id']])
    data['vec'] = [arr for arr in np.concatenate([np.array(vec) for vec in data['vec']])]

    data = pd.DataFrame(data)

    # get families by sizes (largest first)
    sorted_fams = data.groupby([label_name], as_index=False).count()[
        [label_name, 'seq_id']
        ].sort_values(by='seq_id', ascending=False)[label_name].values

    # define labels in binary terms (family vs. not-family)
    for topfam in sorted_fams:
        data[topfam] = [1 if fam == topfam else 0 for fam in data[label_name]]

    return data
