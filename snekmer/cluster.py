"""cluster: Snekmer clustering models.

author: @christinehc

"""
# imports
import numpy as np
from sklearn.base import ClusterMixin
from sklearn.cluster import (
    AgglomerativeClustering,
    DBSCAN,
    Birch,
    OPTICS,
    MiniBatchKMeans,
)
from scipy.cluster.hierarchy import fclusterdata


# wrap scipy correlation clustering into sklearn-like API
class CorrelationClustering(ClusterMixin):
    """Cluster input by correlation using hierarchical clustering.

    Based on `scipy.hierarchy.fcluster`.

    Parameters
    ----------
    **model_params : dict or params
        Parameters for model. See documentation for
        `scipy.hierarchy.fcluster`.

    Attributes
    ----------
    model_params : dict
        Model initialization parameters.
    labels_ : list
        Cluster-assigned labels.

    """

    def __init__(self, **model_params):
        self.model_params = dict(model_params)
        self.labels_ = None

    def fit(self, X):
        self.labels_ = fclusterdata(X, **self.model_params)

    def predict(self, X):
        return fclusterdata(X, **self.model_params)


# define allowed model types
MODELS = {
    "kmeans": MiniBatchKMeans,
    "agglomerative": AgglomerativeClustering,
    "correlation": CorrelationClustering,
    "density": DBSCAN,
    "birch": Birch,
    "optics": OPTICS,
}


# wrapper class for all clustering model types
class KmerClustering(ClusterMixin):
    """Cluster kmer data using various metrics.

    Parameters
    ----------
    method : str
        Clustering method; see `snekmer.cluster.MODELS` for details.
    model_params : type
        Description of parameter `model_params`.

    Attributes
    ----------
    model : type
        Description of attribute `model`.
    method

    """

    def __init__(self, method, model_params={}):
        self.method = method
        self.model = MODELS[method](**model_params)

    def fit(self, X):
        self.model.fit(X)
        self.labels_ = np.full(X.shape[0], -1, dtype=np.intp)

    def predict(self, X):
        return self.model.predict(X)
