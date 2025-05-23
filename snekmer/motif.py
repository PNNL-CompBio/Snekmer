"""motif: Identification of structurally and functionally relevant motifs with Snekmer.
Created on Fri Apr 21 15:25:54 2023

author: @tnitka
"""
# ---------------------------------------------------------
# Imports
# ---------------------------------------------------------
import pandas as pd
import numpy as np


# object to permute training data and retrain
class SnekmerMotif:
    """Permute training data and retrain to find highly distinguishing kmers.

    Parameters
    ----------
    n : int
    Number of permutations to test.
    scores : NDArray
    """

    def __init__(self):
        self.generator = np.random.default_rng()

    def permute(self, X: pd.DataFrame, label, label_col="family"):
        """

        Parameters
        ----------
        X : Dataframe containing matrix of shape (n_kmers, n_features)
            Labeled training data.
        label : str
            Primary family label.
        label_col : str
            Column with family labels.

        Returns
        -------
        Dataframe
            Training data with permuted labels, for retraining and rescoring.

        """
        # save primary family label
        self.primary_label = label
        self.labels = X[label_col].values

        self.generator.shuffle(self.labels)
        X[label_col] = self.labels

        return X

    def p_values(self, X, y: np.ndarray, n: int):
        """

        Parameters
        ----------
        X: Dataframe containing matrix of shape (n_kmers, n_iterations)
            kmer scores from each permutation tested
        y: list or array-like of shape (n_kmers, 1)
            kmer scores from real training data
        n: int
            number of permutations tested

        Returns
        -------
        Dataframe
            matrix containing kmer sequences, scores on real data, number of scores
            on permuted data that exceed that on real data, n_iterations, and
            proportion of scores on permuted data that exceed that on real data.

        """
        self.output_matrix = np.empty((1, 5))
        for i in range(0, len(y) - 1):
            self.seq = X["kmer"].iloc[i]
            self.real_score = y[i]
            self.false_score = X.iloc[i, 1 : (n + 1)].ge(self.real_score).sum()
            self.p = self.false_score / n
            self.vec = np.array(
                [[self.seq, self.real_score, self.false_score, n, self.p]]
            )
            self.output_matrix = np.append(self.output_matrix, self.vec, axis=0)

        else:
            self.output_matrix = np.delete(self.output_matrix, 0, 0)

        self.output = pd.DataFrame(
            self.output_matrix,
            columns=("kmer", "real score", "false positives", "n", "p"),
        )

        return self.output
