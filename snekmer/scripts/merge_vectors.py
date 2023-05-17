# ---------------------------------------------------------
# Imports
# ---------------------------------------------------------
import pickle
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
label = (
    config["score"]["lname"] if str(config["score"]["lname"]) != "None" else "label"
)  # e.g. "family"

# collect all other family tables as negative examples
labels = np.concatenate(
    [
        [skm.utils.split_file_ext(f)[0]]
        * snakemake.params.nseqs[skm.utils.split_file_ext(f)[0]]
        for f in snakemake.input.all_tables
    ],
)
# print(labels)
merged = KmerTable.from_tables(
    [skm.io.load_pickle(f).table for f in snakemake.input.all_tables]
)

# unpack tables using common kmer code set -> combined matrix
X = skm.vectorize.unpack_table(
    merged, sum(snakemake.params.nseqs.values()), kmers=merged.get_kmers()
)

with open(snakemake.output.table, "wb") as f:
    pickle.dump(merged, f)

np.savez_compressed(snakemake.output.labels, labels=labels)
