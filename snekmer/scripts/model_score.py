# ---------------------------------------------------------
# Imports
# ---------------------------------------------------------
import pickle
from datetime import datetime

import pandas as pd
import snekmer as skm
from sklearn.model_selection import train_test_split, StratifiedKFold

# ---------------------------------------------------------
# Files and Parameters
# ---------------------------------------------------------

config = snakemake.config

# ---------------------------------------------------------
# Run script
# ---------------------------------------------------------

# log script start time
start_time = datetime.now()
label = (
    config["score"]["lname"] if str(config["score"]["lname"]) != "None" else "label"
)  # e.g. "family"
with open(snakemake.log[0], "a") as f:
    f.write(f"start time:\t{start_time}\n")

# get kmers for this particular set of sequences
kmers = skm.io.load_pickle(snakemake.input.kmerobj)

# tabulate vectorized seq data
data = list()
kmer_sets = list()
for f in snakemake.input.data:
    kmerlist, df = skm.io.load_npz(f)
    data.append(df)
    kmer_sets.append(kmerlist[0])

# loading all other models and harmonize basis sets
for i in range(len(data)):
    df = data[i]
    kmerlist = kmer_sets[i]
    vecs = skm.utils.to_feature_matrix(df["sequence_vector"].values)
    df["sequence_vector"] = kmers.harmonize(vecs, kmerlist).tolist()

data = pd.concat(data, ignore_index=True)
data["background"] = [f in snakemake.input.bg for f in data["filename"]]

# log conversion step runtime
skm.utils.log_runtime(snakemake.log[0], start_time, step="files_to_df")

# parse family names and only add if some are valid
families = [
    skm.utils.get_family(
        skm.utils.split_file_ext(fn)[0], regex=config["input_file_regex"]
    )
    for fn in data["filename"]
]
if any(families):
    data[label] = families

# binary T/F for classification into family
family = skm.utils.get_family(snakemake.wildcards.nb)
binary_labels = [True if value == family else False for value in data[label]]

# define k-fold split indices
if config["model"]["cv"] > 1:
    cv = StratifiedKFold(n_splits=config["model"]["cv"], shuffle=True)

    # stratify splits by [0,1] family assignment
    for n, (i_train, _) in enumerate(cv.split(data["sequence_vector"], binary_labels)):
        data[f"train_cv-{n + 1:02d}"] = [idx in i_train for idx in data.index]

elif config["model"]["cv"] in [0, 1]:
    i_train, _ = train_test_split(data.index, stratify=binary_labels)
    data["train"] = [idx in i_train for idx in data.index]

# generate family scores and object
scorer = skm.score.KmerScorer()
scorer.fit(
    list(kmers.kmer_set.kmers),
    data,
    skm.utils.get_family(snakemake.wildcards.nb, regex=config["input_file_regex"]),
    label_col=label,
    vec_col="sequence_vector",
    **config["score"]["scaler_kwargs"],
)

# append scored sequences to dataframe
data = data.merge(
    pd.DataFrame(scorer.scores["sample"]), left_index=True, right_index=True
)
if data.empty:
    raise ValueError("Blank df")

# save score loadings
class_probabilities = (
    pd.DataFrame(scorer.probabilities, index=scorer.kmers.basis)
    .reset_index()
    .rename(columns={"index": "kmer"})
)

# log time to compute class probabilities
skm.utils.log_runtime(snakemake.log[0], start_time, step="class_probabilities")

# this will be redundant with the output csv file - except that
#      it will have the harmonized matrix included with it
with open(snakemake.output.matrix, "wb") as f:
    pickle.dump(data, f)

# save all files to respective outputs
delete_cols = ["vec", "sequence_vector"]
for col in delete_cols:
    if col in class_probabilities.columns:
        class_probabilities = class_probabilities.drop(columns=col)
data.drop(columns="sequence_vector").to_csv(
    snakemake.output.data, index=False, compression="gzip"
)

class_probabilities.to_csv(snakemake.output.weights, index=False, compression="gzip")
with open(snakemake.output.scorer, "wb") as f:
    pickle.dump(scorer, f)

# record script endtime
skm.utils.log_runtime(snakemake.log[0], start_time)
