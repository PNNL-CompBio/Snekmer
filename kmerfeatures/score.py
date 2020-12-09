"""score: Scoring and similarity analysis for extracted features.

author: @christinehc
"""
# imports
import pandas as pd

from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics.pairwise import pairwise_distances


# functions
def connection_matrix_from_features(feature_matrix, metric="jaccard"):
    """Calculate similarities between examples based on features
    output from main().

    Parameters
    ----------
    feature_matrix : type
        feature_matrix in the form of rows (proteins), columns (kmers),
        where values are the counts of the kmers for each protein.
    metric : type
        description.

    Returns
    -------
    pandas.DataFrame
        a square matrix with the similarity scores of pairwise
        relationships between proteins.

    """
    # just need to do a pairwise calculation of similarity based on
    #      the specified method. This would be a great place to have
    #      HPC capability because when we get enough examples this is going
    #      to be a lot of calculations
    # I just took this from:
    #   https://pythonpedia.com/en/knowledge-base/37003272/how-to-compute-jaccard-similarity-from-a-pandas-dataframe
    # and I adapted it - not sure if it will work as coded :)

    df = pd.DataFrame(feature_matrix)

    if metric == "jaccard":
        sim = 1 - pairwise_distances(df.T, metric="hamming")
    else:
        sim = pairwise_distances(df.T, metric=metric)

    sim = pd.DataFrame(sim, index=df.columns, columns=df.columns)

    return sim


def cluster_feature_matrix(feature_matrix, method="agglomerative"):
    """Calculate clusters based on the feature matrix.

    Note: sklearn has a wide range of clustering options.

    Parameters
    ----------
    feature_matrix : pandas.DataFrame
        Description of parameter `feature_matrix`.
    method : type
        Description of parameter `method`.

    Returns
    -------
    type
        Description of returned object.

    """
    # another place to have HPC since the clustering can be computationally
    #   intensive. AgglomerativeClustering may not be the best candidate
    #   for parallelization though - not sure

    if method == "agglomerative":
        feat_clusters = AgglomerativeClustering().fit_predict(feature_matrix)

    return feat_clusters
