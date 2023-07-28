# ---------------------------------------------------------
# Imports
# ---------------------------------------------------------
import pickle
from datetime import datetime
from os.path import join

import numpy as np
import pandas as pd
import snekmer as skm
from biotite.sequence.align import KmerTable
from sklearn.model_selection import train_test_split, StratifiedKFold
from sklearn.preprocessing import LabelEncoder

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
kmers = skm.io.load_pickle(snakemake.input.family_table)
family = snakemake.wildcards.nb
# matA = skm.vectorize.unpack_table(tabA, kmerA.n_seqs, kmers=codes)

# collect all other family tables as negative examples
labels = np.load(snakemake.input.labels)["labels"]
labels = np.array([l == "family" for l in labels])
merged = skm.io.load_pickle(snakemake.input.merged_table)
# print(merged)

# for i, f in enumerate(snakemake.input.all_tables, start=1):
#     if snakemake.wildcards.nb not in f:
#         kmerB = skm.io.load_pickle(f)

#         # store values for total set
#         # codes = np.unique(np.concatenate([codes, codesB]))
#         labels = np.concatenate([labels, np.array(["negative"] * kmerB.n_seqs)])
#         table = KmerTable.from_tables([table, kmerB.table])
#         nseqs += kmerB.n_seqs

# unpack tables using common kmer code set -> combined matrix
X = skm.vectorize.unpack_table(
    merged, sum(snakemake.params.nseqs.values()), kmers=merged.get_kmers()
)
print(X.shape, labels.shape)

# generate labels
le = LabelEncoder().fit(labels)
y = le.transform(labels)

# TODO: label bg seqs
# data = pd.concat(data, ignore_index=True)
# data["background"] = [f in snakemake.input.bg for f in data["filename"]]

# log conversion step runtime
skm.utils.log_runtime(snakemake.log[0], start_time, step="combined_matrix")

# binary T/F for classification into family
family = skm.utils.get_family(snakemake.wildcards.nb)
# binary_labels = [True if value == family else False for value in data[label]]

# define k-fold split indices
training_set = dict()
if config["model"]["cv"] > 1:
    cv = StratifiedKFold(n_splits=config["model"]["cv"], shuffle=True)

    # stratify splits by [0,1] family assignment
    for n, (i_train, _) in enumerate(cv.split(X, y)):
        training_set[n] = i_train

elif config["model"]["cv"] in [0, 1]:
    i_train, _ = train_test_split(range(len(y)), stratify=y)
    training_set["train"] = i_train

# generate family scores and object
scorer = skm.score.KmerScorer()
scorer.fit(X, y, family, merged.get_kmers())

# store scored sequences in dataframe
data = scorer.score
if data.empty:
    raise ValueError("Blank df")

# log time to compute class probabilities
skm.utils.log_runtime(snakemake.log[0], start_time, step="probabilities")

# store scored sequences
data.to_csv(snakemake.output.scores, index=False)

# save scorer and kmer probabilityweights
scorer.score.to_csv(snakemake.output.weights, index=False, compression="gzip")
with open(snakemake.output.scorer, "wb") as f:
    pickle.dump(scorer, f)

# record script endtime
skm.utils.log_runtime(snakemake.log[0], start_time)
