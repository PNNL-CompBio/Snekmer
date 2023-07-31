"""model: Machine learning models for class prediction with Snekmer.

author: @christinehc

"""
# imports
import pandas as pd
from typing import Any, Dict, List, Optional
from ._version import __version__
from .vectorize import KmerBasis
from numpy.typing import NDArray
from sklearn.base import BaseEstimator, ClassifierMixin
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier, AdaBoostClassifier
from sklearn.linear_model import LogisticRegression  # LogisticRegressionCV
from sklearn.model_selection import GridSearchCV, cross_validate
from sklearn.pipeline import make_pipeline, Pipeline
from sklearn.svm import SVC

# define default gridsearch param dict
MODEL_PARAMS = {
    #     'scaler__n': [None, 100],
    "clf__class_weight": ["balanced"],
    "clf__min_samples_split": [2, 5, 10, 15, 20, 40],
    "clf__min_samples_leaf": [2, 5, 10, 15, 20],
    "clf__max_depth": [3, 10, 25],
}

NAME2MODEL = {
    "decisiontree": DecisionTreeClassifier,
    "randomforest": RandomForestClassifier,
    "adaboost": AdaBoostClassifier,
    "logistic": LogisticRegression,
    "svc": SVC,
}


# classification models for protein families
class SnekmerModel(BaseEstimator):
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
        self,
        scaler: Optional[Any] = None,
        model: str = "logistic",
        model_params: Optional[Dict[Any, Any]] = {},
        params: Dict[str, List] = MODEL_PARAMS,
        step_name="clf",
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
        # always compute svc probability
        if model == "svc":
            model_params["probability"] = True

        self.scaler = scaler
        self.params = params
        self.model = NAME2MODEL[model](**model_params)
        self.step_name = step_name
        self.snekmer_version = __version__

    def fit(
        self, X: NDArray, y: NDArray, gridsearch: bool = False, verbose: bool = True
    ):
        """Train model using gridsearch-tuned hyperparameters.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            Training input samples.
        y : array-like of shape (n_samples,) or (n_samples, n_outputs)
            Target values (class labels) as integers or strings.
        scores : array-like of shape (n_features,) or (n_features, n_outputs)
            NOT IMPLEMENTED YET-- for KmerScaler integration.
        gridsearch: bool (default: False)
            If True, uses grid search to determine hyperparameters.
        verbose : bool (default: True)
            If True, prints best gridsearch CV results to console.
            Only active if `gridsearch=True`.

        Returns
        -------
        None
            Fits estimator.

        """
        # define pipeline and gridsearch
        if gridsearch:
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
        return self

    def score(self, X: NDArray, y: NDArray):
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

    def predict(self, X: NDArray):
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


class SnekmerModelCV:
    def __init__(self, model: str, cv: int = 5):
        self.model = model
        self.cv = cv

    def train(
        self,
        X: NDArray,
        y: NDArray,
        i: NDArray,
        label: str = "label",
        random_state: Optional[int] = None,
    ):
        """_summary_

        Parameters
        ----------
        X : NDArray
            _description_
        y : NDArray
            _description_
        i : NDArray
            _description_
        label : str, optional
            _description_, by default "label"
        random_state : Optional[int], optional
            _description_, by default None

        Raises
        ------
        ValueError
            _description_

        """
        # set and format input and label arrays; initialize model objects
        cols = [label, "alphabet_name", "k", "scoring"]

        results = {col: [] for col in cols + ["score", "cv_split"]}
        X, y = {i: {} for i in range(self.cv)}, {i: {} for i in range(self.cv)}
        for n in range(self.cv):

            # remove score cols that were generated from full dataset
            unscored_cols = [col for col in list(data.columns) if "_score" not in col]

            # filter data by training data per CV split
            i_train = data[data[f"train_cv-{n + 1:02d}"]].index
            i_test = data[~data[f"train_cv-{n + 1:02d}"]].index
            df_train = data.iloc[i_train][unscored_cols].reset_index(drop=True)
            df_test = data.iloc[i_test][unscored_cols].reset_index(drop=True)
            df_train_labels = [
                True if value == family else False for value in df_train[label]
            ]
            df_test_labels = [
                True if value == family else False for value in df_test[label]
            ]

            # score kmers separately per split
            scorer = skm.score.KmerScorer()
            scorer.fit(
                list(kmer.kmer_set.kmers),
                df_train,
                family,
                label_col=label,
                vec_col="sequence_vector",
                **config["score"]["scaler_kwargs"],
            )

            # append scored sequences to dataframe
            df_train = df_train.merge(
                pd.DataFrame(scorer.scores["sample"]), left_index=True, right_index=True
            )
            if df_train.empty:
                raise ValueError("Blank df")
            df_test = df_test.merge(
                pd.DataFrame(
                    scorer.predict(
                        skm.utils.to_feature_matrix(df_test["sequence_vector"]),
                        list(kmer.kmer_set.kmers),
                    )
                ),
                left_index=True,
                right_index=True,
            ).rename(columns={0: f"{family}_score"})

            # save score loadings
            scores = (
                pd.DataFrame(scorer.probabilities, index=scorer.kmers.basis)
                .reset_index()
                .rename(columns={"index": "kmer"})
            )

            # save X,y array data for plot
            X[n]["train"] = df_train[f"{family}_score"].values.reshape(-1, 1)
            y[n]["train"] = le.transform(df_train_labels).ravel()

            X[n]["test"] = df_test[f"{family}_score"].values.reshape(-1, 1)
            y[n]["test"] = le.transform(df_test_labels).ravel()

        # collate ROC-AUC results
        results[label] += [family] * cv
        results["alphabet_name"] += [alphabet_name.lower()] * cv
        results["k"] += [config["k"]] * cv
        results["scoring"] += ["auc_roc"] * cv
        results["score"] += auc_rocs
        results["cv_split"] += [i + 1 for i in range(cv)]
