"""motif: Identification of structurally and functionally relevant motifs with Snekmer.
Created on Fri Apr 21 15:25:54 2023

author: @tnitka
"""
# ---------------------------------------------------------
# Imports
# ---------------------------------------------------------
import pickle
from datetime import datetime

import snekmer as skm
import pandas as pd
import numpy as np
import snekmer.motif
from typing import Any, Dict, List, Optional
from ._version import __version__
from .vectorize import KmerBasis
from .score import KmerScorer
from .model import SnekmerModel, SnekmerModelCV
#from numpy.typing import NDArray
from sklearn.base import BaseEstimator, ClassifierMixin
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier, AdaBoostClassifier
from sklearn.linear_model import LogisticRegression  # LogisticRegressionCV
from sklearn.model_selection import GridSearchCV, cross_validate
from sklearn.pipeline import make_pipeline, Pipeline
from sklearn.svm import SVC

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
        self.scorer = skm.score.KmerScorer()
    
    def permute(self, X: pd.DataFrame, label, label_col="family"):
        """
        
        Parameters
        ----------
        X : Dataframe containing matrix of shape (n_samples, n_features)
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
        
        #self.permuted_data = X
        self.permuted_labels = self.generator.permutation(self.labels)
        #self.mask = np.zeros_like(X, dtype=bool)
        #self.mask[:, 1] = True
        #np.place(self.permuted_data, self.mask, self.permuted_labels)
        self.permuted_data = X.assign(label_col=self.permuted_labels)
        
        return self.permuted_data
        
    def p_values(self, X, y: np.ndarray, n: int):
        """
        
        Parameters
        ----------
        X: Dataframe containing matrix of shape (n_samples, n_iterations)
            kmer scores from each permutation tested
        y: list or array-like of shape (n_samples, 1)
            kmer scores from real training data
        n: int
            number of permutations tested

        Returns
        -------
        NDArray
            matrix containing kmer sequences, scores on real data, number of scores
            on permuted data that exceed that on real data, n_iterations, and
            proportion of scores on permuted data that exceed that on real data.

        """
        
        self.output = pd.DataFrame(columns=('kmer', 'real score', 'false positives', 'n', 'p'))
        for i in range(1, len(y)):
            self.seq = self.labels[i]
            self.real_score = y[i]
            # self.score_list = X.iloc[i, :].values.tolist()
            # self.false_score = sum(j > self.real_score for j in pd.to_numeric(X.iloc[i]))
            self.false_score = X.iloc[i, 1:len(y)].gt(self.real_score).sum()
            self.p = self.false_score/n
            self.dict = {
                 'kmer': [self.seq],
                 'real score': [self.real_score],
                 'false positives': [self.false_score],
                 'n': [n],
                 'p': [self.p]}
            self.vec = pd.DataFrame.from_dict(self.dict)
            self.output.merge(self.vec, on='kmer')

            
        # else:
        #     self.output = np.delete(self.output, 1, 0)
        
        return self.output