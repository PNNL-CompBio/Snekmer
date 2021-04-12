"""model: Machine learning models for Kmer family prediction.

author: @christinehc
"""
# imports
import json
import numpy as np
import pandas as pd
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier, AdaBoostClassifier
from sklearn.model_selection import GridSearchCV
from sklearn.pipeline import make_pipeline, Pipeline
from .utils import get_family

# define default gridsearch param dict
MODEL_PARAMS = {
    'clf__class_weight': [None, 'balanced'],
    'clf__min_samples_split': [2, 5, 10, 15, 20, 40],
    'clf__min_samples_leaf': [2, 5, 10, 15, 20],
    'clf__max_depth': [3, 10, 25],
  }


# classes
class KmerScaler:
    """Scale and reduce kmer vector based on scores.

    Parameters
    ----------
    model : type
        Description of parameter `model`.

    Attributes
    ----------
    cv : type
        Description of attribute `cv`.
    model_params : type
        Description of attribute `model_params`.

    """

    def __init__(self, scores, n=100):
        """Initialize KmerScaler object.

        Attributes
        ----------
        kmer_basis_idx : numpy.ndarray
            Description of parameter `kmer_basis_idx`.
        kmer_basis_score : dict
            Description of parameter `kmer_basis_score`.

        Returns
        -------
        type
            Description of returned object.

        """
        self.n = n
        self.input_shape = None
        self.kmer_basis_idx = None
        self.kmer_basis_score = dict()
        self.scores = scores


    def fit(self, X, threshold=None):
        """Fit scaler based on list of scores.

        Parameters
        ----------
        scores : type
            Description of parameter `scores`.
        threshold : type
            Description of parameter `threshold`.
        n : int

        Returns
        -------
        type
            Description of returned object.

        """
        if not isinstance(scores, np.ndarray):
            scores = np.array(scores)
        self.input_shape = scores.shape
        if threshold is not None:
            indices = np.where(np.array(scores > threshold))[0]
        elif self.n is not None:
            indices = scores.ravel().argsort()[:-self.n - 1:-1]
        else:
            raise ValueError("One of either `threshold` or `n` must"
                             " be specified.")
        self.kmer_basis_idx = indices
        self.kmer_basis_score = [(i, score) for i, score in zip(
            self.kmer_basis_idx, scores[self.kmer_basis_idx])]
        return

    def transform(self, array):
        """Reduce vector to kmer basis set.

        Parameters
        ----------
        array : type
            Description of parameter `array`.

        Returns
        -------
        type
            Description of returned object.

        """
        if self.kmer_basis_idx is None or not any(self.kmer_basis_idx):
            raise AttributeError("Kmer basis set not defined;"
                                 " must fit scores to scaler.")
        if not isinstance(array, np.ndarray):
            array = np.array(array)
        if array.shape != self.input_shape:
            if array[0].shape != self.input_shape:
                raise ValueError(f"Input vector shape {array.shape}"
                                 f" does not match fitted vector shape"
                                 f" {self.input_shape}.")
            return np.array([a[self.kmer_basis_idx] for a in array])
        return array[self.kmer_basis_idx]


# define default gridsearch param dict
MODEL_PARAMS = {
#     'scaler__n': [None, 100],
    'clf__class_weight': [None, 'balanced'],
    'clf__min_samples_split': [2, 5, 10, 15, 20, 40],
    'clf__min_samples_leaf': [2, 5, 10, 15, 20],
    'clf__max_depth': [3, 10, 25],
  }


class KmerModel:
    """Classify a protein family using kmer vectors as input.

    Attributes
    ----------
    input_shape : numpy.ndarray
        Shape of input kmer feature vector.
    kmer_basis_idx : type
        Description of attribute `kmer_basis_idx`.
    kmer_basis_score : dict
        Description of attribute `kmer_basis_score`.

    """
    def __init__(self,
                 scaler=None,
                 model=DecisionTreeClassifier(),
                 params=MODEL_PARAMS,
                 step_name="clf"):
        """Initialize KmerModel object.

        Parameters
        ----------
        model : type
            Description of parameter `model`.
        params : type
            Description of parameter `params`.
        step_name : type
            Description of parameter `step_name`.

        Returns
        -------
        type
            Description of returned object.

        """
        self.scaler = scaler
        self.params = params
        self.model = model
        self.step_name = step_name

        self.pipeline = None
        self.search = None

    def fit(self, x_train, y_train, verbose=True):
        """Train model using gridsearch-tuned hyperparameters.

        Parameters
        ----------
        x_train : type
            Description of parameter `x_train`.
        y_train : type
            Description of parameter `y_train`.
        scores : type
            Description of parameter `scores`.
        verbose : type
            Description of parameter `verbose`.

        Returns
        -------
        type
            Description of returned object.

        """
        # define pipeline and gridsearch
        self.pipeline = Pipeline(steps=[("scaler", self.scaler), (self.step_name, self.model)])
        self.search = GridSearchCV(self.pipeline, self.params)

        # use gridsearch on training data to find best parameter set
        self.search.fit(x_train, y_train)
        if verbose:
            print("best mean cross-validation score: {:.3f}".format(
                self.search.best_score_
                ))
            print("best parameters:", self.search.best_params_)

        # fit model with best gridsearch parameters
#         self.scaler.fit(scores)
        self.model.fit(x_train, y_train)
        return

    def score(self, x_test, y_test):
        """Score model on test set.

        Parameters
        ----------
        x : type
            Description of parameter `x`.
        y : type
            Description of parameter `y`.

        Returns
        -------
        type
            Description of returned object.

        """
        return self.search.score(x_test, y_test)


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


# def train_model(model=DecisionTreeClassifier, n_splits=5):
#
#     # define k-fold cross validation
#     skfold = StratifiedKFold(n_splits=n_splits, shuffle=True)
#     scoring = ['accuracy', 'precision_macro',
#                'recall_macro', 'f1_macro', 'roc_auc']
#     score_dict_keys = {f"test_{score_method}": [] for score_method in scoring}
#
#     cv_results = {'family': [], 'norm': [],
#                   'model_params': [], 'model': [],
#                   # 'X_train': [],
#                   # 'X_test': [],
#                   # 'y_train': [],
#                   # 'y_test': [],
#                   **score_dict_keys}
#     for fam in sorted_fams[:10]:
#         Xf, yf, Yf = np.asarray([arr for arr in data['vec'].values]), data[fam].values, data['family'].values
#         Xf = (Xf > 0) * 1
#         # Xfs = StandardScaler().fit_transform(Xf)
#         Xf_train, Xf_test, yf_train, yf_test = train_test_split(
#             Xf, yf, test_size=.2, random_state=9, stratify=yf
#         )
#
#         param_grid = {'n_estimators': [10, 50, 100, 200],
#                       'criterion': ['gini', 'entropy'],
#                       'max_depth': [3, 5, 10, 15],
#                       'min_samples_split': [2, 5, 10, 15, 20, 25, 30, 35, 40],
#                       'min_samples_leaf': [2, 5, 10, 15, 20],
#                       'max_features': [None, 'auto', 'sqrt', 'log2'],
#                       'class_weight': [None, 'balanced']}
#         gridb = GridSearchCV(RandomForestClassifier(),
#                              param_grid=param_grid)
#         gridb.fit(Xf_train, yf_train)
#         print("best mean cross-validation score: {:.3f}".format(gridb.best_score_))
#         print("best parameters:", gridb.best_params_)
#         print("test-set score: {:.3f}".format(gridb.score(Xf_test, yf_test)))
#         # best_params_b = {'criterion': 'entropy',
#         #                  'max_depth': 5,
#         #                  'min_samples_split': 2,
#         #                  'min_samples_leaf': 2,
#         #                  'max_features': None,
#         #                  'class_weight': 'balanced'}  # gridb.best_params_
#         best_params_b = gridb.best_params_
#         return
