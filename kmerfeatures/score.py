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


# For Christine
#   This code will calculate the probabilities for features being
#       in a class (as defined by example_classes)
# Input: feature matrix and an example_classes vector which has as
#        many entries as the number of examples and indicate a class
def feature_class_probabilities(feature_matrix, example_labels):
    # coding this first for the binary case - in class or not

    # check to make sure the labels are the same size as number of examples
    # check labels to make sure they contain two classes
    #  TBD

    # convert the feature count matrix into binary presence/absence
    feature_matrix = (feature_matrix>0)*1

    colnames = feature_matrix.keys()

    #results = pd.DataFrame(data=colnames, columns=["CountPos","ProbPos","CountNeg","ProbNeg","Score"])
    results = pd.DataFrame(index=colnames, columns=["CountPos","ProbPos","CountNeg","ProbNeg","Score"])

    pos_tot = sum(example_labels)
    neg_tot = len(example_labels) - pos_tot

    # for every feature in the input matrix
    for key in feature_matrix.keys():
        features = feature_matrix[key]
        pos_score = example_labels * features
        neg_score = ((example_labels==0)*1) * features

        # probability that the presence of this kmer maps to the positive
        #   examples
        pos_score_p = sum(pos_score)/pos_tot
        neg_score_p = sum(neg_score)/neg_tot

        results.loc[key] = [sum(pos_score), pos_score_p, sum(neg_score), neg_score_p, pos_score_p-neg_score_p]

    return(results)
