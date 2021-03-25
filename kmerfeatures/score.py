"""score: Scoring and similarity analysis for extracted features.

author: @christinehc
"""
# imports
import pandas as pd

from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics.pairwise import pairwise_distances


# functions
def connection_matrix_from_features(feature_matrix, metric="jaccard"):
    """Calculate similarities based on features output from main().

    Parameters
    ----------
    feature_matrix : numpy.ndarray
        feature_matrix in the form of rows (proteins), columns (kmers),
        where values are the counts of the kmers for each protein.
    metric : str
        description.

    Returns
    -------
    pandas.DataFrame
        a square matrix with the similarity scores of pairwise
        relationships between proteins.

    """
    if metric == "jaccard":
        sim = 1 - pairwise_distances(feature_matrix.T, metric="hamming")
    else:
        sim = pairwise_distances(feature_matrix.T, metric=metric)
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
        clusters = AgglomerativeClustering().fit_predict(feature_matrix)

    return clusters


def feature_class_probabilities(feature_matrix, example_labels):
    """Calculate probabilities for features being in a defined class.

    Parameters
    ----------
    feature_matrix : type
        Feature matrix.
    example_labels : type
        Example classes vector which has as many entries as the
        number of examples and indicates a class.

    Returns
    -------
    type
        Description of returned object.

    """
    # coding this first for the binary case - in class or not

    # check to make sure the labels are the same size as number of examples
    # check labels to make sure they contain two classes
    #  TBD

    # convert the feature count matrix into binary presence/absence
    feature_matrix = (feature_matrix > 0) * 1

    cols = feature_matrix.keys()

    # results = pd.DataFrame(data=colnames, columns=["CountPos","ProbPos","CountNeg","ProbNeg","Score"])
    results = pd.DataFrame(index=cols,
                           columns=["CountPos", "ProbPos",
                                    "CountNeg", "ProbNeg", "Score"])

    pos_total = sum(example_labels)
    neg_total = len(example_labels) - pos_total

    # for every feature in the input matrix
    for key in feature_matrix.keys():
        features = feature_matrix[key]
        pos_score = example_labels * features
        neg_score = ((example_labels == 0) * 1) * features

        # probability that presence of this kmer maps to positive examples
        pos_prob_score = sum(pos_score) / pos_total
        neg_prob_score = sum(neg_score) / neg_total

        results[key] = [sum(pos_score), pos_prob_score,
                        sum(neg_score), neg_prob_score,
                        pos_prob_score - neg_prob_score]

    return results
