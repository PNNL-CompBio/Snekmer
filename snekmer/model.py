"""model: Machine learning models for class prediction with Snekmer.

author: @christinehc

"""
# imports
import json

import numpy as np
import pandas as pd

from .utils import get_family, to_feature_matrix
from .score import (
    feature_class_probabilities,
    apply_feature_probabilities,
    KmerScoreScaler,
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
    "clf__class_weight": ["balanced"],
    "clf__min_samples_split": [2, 5, 10, 15, 20, 40],
    "clf__min_samples_leaf": [2, 5, 10, 15, 20],
    "clf__max_depth": [3, 10, 25],
}

MODEL_NAME = {
    "decisiontree": DecisionTreeClassifier(),
    "randomforest": RandomForestClassifier(),
    "adaboost": AdaBoostClassifier(),
}


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

    def __init__(
        self, scaler=None, model="decisiontree", params=MODEL_PARAMS, step_name="clf"
    ):
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
        self.pipeline = Pipeline(
            steps=[("scaler", self.scaler), (self.step_name, self.model)]
        )
        self.search = GridSearchCV(self.pipeline, self.params)

        # use gridsearch on training data to find best parameter set
        self.search.fit(X, y)
        if verbose:
            print(
                "best mean cross-validation score: {:.3f}".format(
                    self.search.best_score_
                )
            )
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
def format_data_df(filenames, label_name="family"):
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
    data = {"seq_id": [], "vec": [], label_name: []}
    for fn in filenames:
        with open(fn, "r") as f:
            tmp = json.load(f)
            data["seq_id"] += [tmp["seq_id"]]
            data["vec"] += [np.array(tmp["vector"])]
            data[label_name] += [get_family(fn)]

    data[label_name] = np.hstack(
        np.array(
            [[label] * len(data["vec"][i]) for i, label in enumerate(data[label_name])],
            dtype="object",
        )
    )
    data["seq_id"] = np.hstack([np.array(s) for s in data["seq_id"]])
    data["vec"] = [
        arr for arr in np.concatenate([np.array(vec) for vec in data["vec"]])
    ]

    data = pd.DataFrame(data)

    # get families by sizes (largest first)
    sorted_fams = (
        data.groupby([label_name], as_index=False)
        .count()[[label_name, "seq_id"]]
        .sort_values(by="seq_id", ascending=False)[label_name]
        .values
    )

    # define labels in binary terms (family vs. not-family)
    for topfam in sorted_fams:
        data[topfam] = [1 if fam == topfam else 0 for fam in data[label_name]]

    return data
