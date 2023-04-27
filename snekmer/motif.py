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
    
    def __init__(self, n=2000):
        self.scores = skm.io.load_pickle(snakemake.input.weights)
        self.real_data = skm.io.load_pickle(snakemake.input.matrix)
        self.labels = skm.io.load_pickle(snakemake.input.data)
        self.generator = np.random.default_rng()
        self.scorer = skm.score.KmerScorer()
    
    def permute(self, X: np.ndarray, y: np.ndarray):
        """
        
        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            Labeled training data.
        y : list or array-like of shape (n_samples, 1)
            Family labels from training data.

        Returns
        -------
        NDArray
            Training data with permuted labels, for retraining and rescoring.

        """
        
        self.permuted_labels = self.generator.permutation(y)
        self.unlabeled = np.delete(X, 0, 1)
        self.permuted_data = np.hstack(self.permuted_labels, self.unlabeled)
        
        return self.permuted_data
        
    def p_values(self, X: np.ndarray, y: np.ndarray, n: int):
        """
        
        Parameters
        ----------
        X: array-like of shape (n_samples, n_iterations)
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
        
        self.output = np.empty((1,5))
        for i in range(len(y)):
            self.seq = self.labels[i]
            self.real_score = self.scores[i]
            self.false_score = sum(j > self.real_score for j in X[i, :])
            self.p = self.false_score/n
            self.vec = [self.seq, self.real_score, self.false_score, n, self.p],
            self.output = np.append(self.output, self.vec, 0)
            
        else:
            self.output = np.delete(self.output, 1, 0)
        
        return self.output