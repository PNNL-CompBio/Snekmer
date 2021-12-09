"""plot: Plotting functions for kmer pipeline

author: @christinehc
"""

# https://stackoverflow.com/questions/29656550/how-to-plot-pr-curve-over-10-folds-of-cross-validation-in-scikit-learn
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import (accuracy_score,
                             auc,
                             average_precision_score,
                             roc_curve,
                             precision_recall_curve)
from sklearn.model_selection import (KFold,
                                     train_test_split,
                                     RandomizedSearchCV,
                                     StratifiedKFold)
from sklearn.svm import SVC
from sklearn.preprocessing import LabelEncoder, StandardScaler
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE


def show_cv_roc_curve(clf, cv, X, y, title='ROC Curve', ax=None, dpi=400):
    """Plot cross-validated receiver operator characteristics curve.

    Adapted from example in sklearn documentation [1].
    [1]: http://scikit-learn.org/stable/auto_examples/model_selection/plot_roc_crossval.html#sphx-glr-auto-examples-model-selection-plot-roc-crossval-py

    Parameters
    ----------
    clf : Classifer object
        Classifier object (e.g. DecisionTreeClassifier()).
    cv : sklearn.model_selection.StratifiedKFold object
        Cross-validation object.
    X : pandas.DataFrame
        Feature dataframe.
    y : pandas.DataFrame
        Response dataframe.
    title : str
        Plot title.
    ax : matplotlib.axes.Axes object or None (default: None)
        Axis object, if already created.
    dpi : int (default: 400)
        Figure resolution.

    Returns
    -------
    None

    """
    # initialize figure and define all axes
    fig = plt.figure(dpi=dpi, figsize=(8, 4), facecolor='white')
    ax = ax or plt.gca()

    # create empty arrays to store data
    tprs, aucs = [], []
    mean_fpr = np.linspace(0, 1, 100)

    # take each cv result
    for i, (train, test) in enumerate(cv.split(X, y)):
        probabilities = clf.fit(X[train], y[train]).predict_proba(X[test])
        # calculate roc curve/area, and interpolate values
        fpr, tpr, _ = roc_curve(y[test], probabilities[:, 1])
        tprs.append(np.interp(mean_fpr, fpr, tpr))

        tprs[-1][0] = 0.0
        roc_auc = auc(fpr, tpr)
        aucs.append(roc_auc)
        ax.plot(fpr, tpr, lw=1, alpha=0.3,
                label=f"ROC fold {i} (AUC = {roc_auc:0.2f})")

    # plot y = x (50% chance) reference line
    ax.plot([0, 1], [0, 1], linestyle='--', lw=2, color='r',
            label='Chance', alpha=0.8)

    # calculate mean and stdev roc_auc and show in plot
    mean_tpr = np.mean(tprs, axis=0)
    mean_tpr[-1] = 1.0
    mean_auc = auc(mean_fpr, mean_tpr)
    std_auc = np.std(aucs)
    ax.plot(mean_fpr, mean_tpr, color='b',
            label=(f"Mean ROC (AUC = {mean_auc:0.2f}"
                   r" $\pm$ "
                   f" {std_auc:0.2f})"),
            lw=2, alpha=.8)

    # draw shaded area for +- 1 stdev from mean auc roc
    std_tpr = np.std(tprs, axis=0)
    tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
    tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
    ax.fill_between(mean_fpr, tprs_lower, tprs_upper, color='grey', alpha=.2,
                    label=r'$\pm$ 1 std. dev.')

    # set fixed plot limits from ~0 to ~1
    ax.set_xlim([-0.05, 1.05])
    ax.set_ylim([-0.05, 1.05])
    ax.set_xlabel('False Positive Rate')
    ax.set_ylabel('True Positive Rate')
    ax.set_title(title)
    ax.legend(bbox_to_anchor=(1.05, 0.5), loc='center left', borderaxespad=0.)

    return fig, ax


def show_cv_pr_curve(clf, cv, X, y, title='PR Curve', ax=None, dpi=400):
    """Plot cross-validated precision-recall curve.

    Parameters
    ----------
    clf : Classifer object
        Classifier object (e.g. DecisionTreeClassifier()).
    cv : sklearn.model_selection.StratifiedKFold object
        Cross-validation object.
    X : pandas.DataFrame
        Feature dataframe.
    y : pandas.DataFrame
        Response dataframe.
    title : str
        Plot title.
    ax : matplotlib.axes.Axes object or None (default: None)
        Axis object, if already created.
    dpi : int (default: 400)
        Figure resolution.

    Returns
    -------
    None

    """
    # initialize figure and define all axes
    fig = plt.figure(dpi=dpi, figsize=(8, 4), facecolor='white')
    ax = ax or plt.gca()

    # create empty arrays to store data
    y_real, y_proba = [], []

    for i, (train, test) in enumerate(cv.split(X, y)):
        probabilities = clf.fit(X[train], y[train]).predict_proba(X[test])
        # Compute ROC curve and area the curve
        precision, recall, _ = precision_recall_curve(y[test], probabilities[:, 1])
        avg_precision = average_precision_score(y[test], probabilities[:, 1])

        # Plotting each individual PR Curve
        ax.plot(recall, precision, lw=1, alpha=0.3,
                label=f"PR fold {i} (AUC = {avg_precision:0.2f})")

        y_real.append(y[test])
        y_proba.append(probabilities[:, 1])

    y_real = np.concatenate(y_real)
    y_proba = np.concatenate(y_proba)

    precision, recall, _ = precision_recall_curve(y_real, y_proba)
    avg_total_precision = average_precision_score(y_real, y_proba)

    # plot average p-r curve
    ax.plot(recall, precision, color='b',
            label=f"Precision-Recall (AUC = {avg_total_precision:0.2f})",
            lw=2, alpha=.8)

    # set fixed plot limits from ~0 to ~1
    ax.set_xlim([-0.05, 1.05])
    ax.set_ylim([-0.05, 1.05])
    ax.set_xlabel('Recall')
    ax.set_ylabel('Precision')
    ax.set_title(title)
    ax.legend(bbox_to_anchor=(1.05, 0.5), loc='center left', borderaxespad=0.)
    return fig, ax


# plot PCA explained variance
def show_explained_variance_curve(feature_matrix):
    X_scaled = StandardScaler().fit_transform(feature_matrix)
    pca = PCA()
    pca.fit(X_scaled)

    fig, ax = plt.subplots(dpi=200, figsize=(4, 3))
    ax.plot(
        [0] + list(np.arange(1, pca.n_components_ + 1)),
        [0] + list(np.cumsum(pca.explained_variance_ratio_))
    )
    ax.set_ylabel("PCA explained variance ratio")

    return fig, ax


# plot clusters using tsne (requires scaling + pca)
def get_tsne_clusters(feature_matrix, clusters, **tsne_args):
    X_scaled = StandardScaler().fit_transform(feature_matrix)
    X_pca = PCA().fit_transform(X_scaled)
    X_embedded = TSNE(**tsne_args).fit_transform(X_pca)#[:, :100])#_unscaled)

    fig, ax = plt.subplots(dpi=200, figsize=(4, 3))
    sns.scatterplot(x=X_embedded[:, 0], y=X_embedded[:, 1], hue=clusters, ax=ax)
    return fig, ax
