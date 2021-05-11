"""score: Scoring and similarity analysis for extracted features.

author: @christinehc / @biodataganache
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
    numpy.ndarray
        a square matrix with the similarity scores of pairwise
        relationships between proteins.

    """
    if metric == "jaccard":
        sim = 1 - pairwise_distances(feature_matrix, metric="hamming")  # .T?
    else:
        sim = pairwise_distances(feature_matrix, metric=metric)
    return sim


def cluster_feature_matrix(feature_matrix, method="agglomerative", **kwargs):
    """Calculate clusters based on the feature matrix.

    Note: sklearn has a wide range of clustering options.

    Parameters
    ----------
    feature_matrix : pandas.DataFrame
        Description of parameter `feature_matrix`.
    method : type
        Description of parameter `method`.
    **kwargs : dict
        Keyword arguments for clustering class.

    Returns
    -------
    type
        Description of returned object.

    """
    # another place to have HPC since the clustering can be computationally
    #   intensive. AgglomerativeClustering may not be the best candidate
    #   for parallelization though - not sure

    if method == "agglomerative":
        clusters = AgglomerativeClustering(n_clusters=n_clusters).fit_predict(feature_matrix)

    return clusters


def feature_class_probabilities(feature_matrix, example_labels, df=True):
    """Calculate probabilities for features being in a defined class.

    Parameters
    ----------
    feature_matrix : type
        Feature matrix.
    example_labels : type
        Example classes vector which has as many entries as the
        number of examples and indicates a class.
    df : bool
        If True, returns output as a pandas DataFrame;
        if False, returns output as a dictionary (default: True).

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
    feature_matrix = (feature_matrix > 0) * 1.0

    cols = ["count_pos", "prob_pos", "count_neg", "prob_neg", "score"]

    # count label totals in each category
    pos_total = sum(example_labels)
    neg_total = len(example_labels) - pos_total

    results = {key: [] for key in cols}
    # for every feature in the input matrix
    for i, features in enumerate(feature_matrix):
        pos_score = example_labels[i] * features
        neg_score = ((example_labels[i] == 0) * 1) * features

        # probability that presence of this kmer maps to positive examples
        pos_prob_score = sum(pos_score) / pos_total
        neg_prob_score = sum(neg_score) / neg_total

        results['count_pos'] = results['count_pos'] + [sum(pos_score)]
        results['prob_pos'] = results['prob_pos'] + [pos_prob_score]
        results['count_neg'] = results['count_neg'] + [sum(neg_score)]
        results['prob_neg'] = results['prob_neg'] + [neg_prob_score]
        results['score'] = results['score'] + [pos_prob_score - neg_prob_score]

    if df:
        return pd.DataFrame(results)

    return results

# This code will apply the feature class probabilities derived in the previous functions
#   to a new set of examples as a weighted score. The idea here is to use a very simple
#   approach to classification.
def apply_feature_probabilities(feature matrix, probabilities):
    """Calculate probabilities for features being in a defined class.

    Parameters
    ----------
    feature_matrix : type
        Feature matrix.
    probabilities : type
        Pandas data frame of kmer probabilities from "feature_class_probabilities".

    Returns
    -------
    type
        Pandas data frame of examples (proteins) with associated score

    """


    # first make sure the feature matrix is binary

    # TODO: we need to check that the probabilities are the same
    #       as the input feature matrix - that is, that they
    #       refer to the same features in the same order

    feature_matrix = (feature_matrix>0)*1.0
    scores = {'score'}

    for i, features in enumerate(feature_matrix):
        # for a first attempt we can just multiply the probability
        #    value by the kmer presence vector which will give us
        #    a weighted score. No idea if this will work well or not
        score = features * probabilities['score']
        scores[i] = score

    return(pd.DataFrame(scores))
