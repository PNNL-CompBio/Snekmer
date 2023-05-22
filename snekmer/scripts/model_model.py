# ---------------------------------------------------------
# Imports
# ---------------------------------------------------------

import pickle
from datetime import datetime
from os import makedirs
from os.path import exists, join

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import snekmer as skm
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split, StratifiedKFold
from sklearn.preprocessing import LabelEncoder

# ---------------------------------------------------------
# Files and Parameters
# ---------------------------------------------------------

config = snakemake.config

# change matplotlib backend to non-interactive
plt.switch_backend("Agg")

# ---------------------------------------------------------
# Run script
# ---------------------------------------------------------

# create lookup table to match vector to sequence by file+ID
# lookup = {}
# for f in input.raw:
#    loaded = np.load(f)
#    lookup.update(
#        {
#            (splitext(basename(f))[0], seq_id): seq_vec
#            for seq_id, seq_vec in zip(loaded["ids"], loaded["vecs"])
#        }
#    )
# print(lookup)
with open(snakemake.input.matrix, "rb") as f:
    data = pickle.load(f)

# load all input data and encode rule-wide variables
# data = pd.read_csv(input.data)
# data["sequence_vector"] = [
#    lookup[(seq_f, seq_id)]
#    for seq_f, seq_id in zip(data["filename"], data["sequence_id"])
# ]

scores = pd.read_csv(snakemake.input.weights)
family = skm.utils.get_family(
    skm.utils.split_file_ext(snakemake.input.weights)[0],
    regex=config["input_file_regex"],
)
# get kmers for this particular set of sequences
with open(snakemake.input.kmerobj, "rb") as f:
    kmer = pickle.load(f)

cv = config["model"]["cv"]

# set category label name (e.g. "family")
label = config["score"]["lname"] if str(config["score"]["lname"]) != "None" else "label"

# prevent kmer NA being read as np.nan
if config["k"] == 2:
    scores["kmer"] = scores["kmer"].fillna("NA")

# get alphabet name
if config["alphabet"] in skm.alphabet.ALPHABET_ORDER.keys():
    alphabet_name = skm.alphabet.ALPHABET_ORDER[config["alphabet"]].capitalize()
else:
    alphabet_name = str(config["alphabet"]).capitalize()

# generate [0, 1] labels for binary family assignment
binary_labels = [True if value == family else False for value in data[label]]
le = LabelEncoder()
le.fit(binary_labels)

# set and format input and label arrays; initialize model objs
# need to make training or test set?
X_all = data[f"{family}_score"].values.reshape(-1, 1)
y_all = le.transform(binary_labels).ravel()

# set random seed if specified
rng = np.random.default_rng()
random_state = (
    config["model"]["random_state"]
    if str(config["model"]["random_state"]) != "None"
    else rng.integers(low=0, high=32767)  # max for int16
)

# set and format input and label arrays; initialize model objs
cols = [label, "alphabet_name", "k", "scoring"]
results = {col: [] for col in cols + ["score", "cv_split"]}
X, y = {i: {} for i in range(cv)}, {i: {} for i in range(cv)}
for n in range(cv):

    # remove score cols that were generated from full dataset
    unscored_cols = [col for col in list(data.columns) if "_score" not in col]

    # filter data by training data per CV split
    i_train = data[data[f"train_cv-{n + 1:02d}"]].index
    i_test = data[~data[f"train_cv-{n + 1:02d}"]].index
    df_train = data.iloc[i_train][unscored_cols].reset_index(drop=True)
    df_test = data.iloc[i_test][unscored_cols].reset_index(drop=True)
    df_train_labels = [True if value == family else False for value in df_train[label]]
    df_test_labels = [True if value == family else False for value in df_test[label]]

    # score kmers separately per split
    scorer = skm.score.KmerScorer()
    scorer.fit(
        list(kmer.kmer_set.kmers),
        df_train,
        family,
        label_col=label,
        vec_col="sequence_vector",
        **config["score"]["scaler_kwargs"],
    )

    # append scored sequences to dataframe
    df_train = df_train.merge(
        pd.DataFrame(scorer.scores["sample"]), left_index=True, right_index=True
    )
    if df_train.empty:
        raise ValueError("Blank df")
    df_test = df_test.merge(
        pd.DataFrame(
            scorer.predict(
                skm.utils.to_feature_matrix(df_test["sequence_vector"]),
                list(kmer.kmer_set.kmers),
            )
        ),
        left_index=True,
        right_index=True,
    ).rename(columns={0: f"{family}_score"})

    # save score loadings
    scores = (
        pd.DataFrame(scorer.probabilities, index=scorer.kmers.basis)
        .reset_index()
        .rename(columns={"index": "kmer"})
    )

    # save X,y array data for plot
    X[n]["train"] = df_train[f"{family}_score"].values.reshape(-1, 1)
    y[n]["train"] = le.transform(df_train_labels).ravel()

    X[n]["test"] = df_test[f"{family}_score"].values.reshape(-1, 1)
    y[n]["test"] = le.transform(df_test_labels).ravel()

# ROC-AUC figure
clf = skm.model.SnekmerModel(
    model="logistic",
    params={
        "random_state": random_state,
        "solver": "liblinear",
        "class_weight": "balanced",
    },
)
# clf = LogisticRegression(
#     random_state=random_state, solver="liblinear", class_weight="balanced"
# )
fig, ax, auc_rocs = skm.plot.cv_roc_curve(
    clf, X, y, title=f"{family} ROC Curve ({alphabet_name}, k={config['k']})"
)

# collate ROC-AUC results
results[label] += [family] * cv
results["alphabet_name"] += [alphabet_name.lower()] * cv
results["k"] += [config["k"]] * cv
results["scoring"] += ["auc_roc"] * cv
results["score"] += auc_rocs
results["cv_split"] += [i + 1 for i in range(cv)]

# save ROC-AUC figure
if not exists(snakemake.output.figs):
    makedirs(snakemake.output.figs)
fig.savefig(
    join(
        snakemake.output.figs,
        (f"{family}_roc-auc-curve_{alphabet_name.lower()}" f"_k-{config['k']:02d}.png"),
    ),
    dpi=fig.dpi,
)
fig.clf()
plt.close("all")

# PR-AUC figure
fig, ax, pr_aucs = skm.plot.cv_pr_curve(
    clf, X, y, title=f"{family} PR Curve ({alphabet_name}, k={config['k']})"
)

# collate PR-AUC results
results[label] += [family] * cv
results["alphabet_name"] += [alphabet_name.lower()] * cv
results["k"] += [config["k"]] * cv
results["scoring"] += ["pr_auc"] * cv
results["score"] += pr_aucs
results["cv_split"] += [i + 1 for i in range(cv)]

# save PR-AUC figure
fig.savefig(
    join(
        snakemake.output.figs,
        (f"{family}_aupr-curve_{alphabet_name.lower()}" f"_k-{config['k']:02d}.png"),
    ),
    dpi=fig.dpi,
)
fig.clf()
plt.close("all")

# save model
clf.fit(X_all, y_all)
with open(snakemake.output.model, "wb") as save_model:
    pickle.dump(clf, save_model)

# save full results
pd.DataFrame(results).to_csv(snakemake.output.results, index=False)
