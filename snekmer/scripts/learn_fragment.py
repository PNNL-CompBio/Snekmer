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

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.csv as csv
import seaborn as sns
import sklearn
from Bio import SeqIO
from scipy.interpolate import interp1d
from scipy.stats import rankdata
import random
import snekmer as skm
from scipy.ndimage import gaussian_filter1d
import sklearn.metrics.pairwise
import itertools
import os
import shutil
import sys
# ---------------------------------------------------------
# Files and Parameters
# ---------------------------------------------------------

config = snakemake.config

# change matplotlib backend to non-interactive
plt.switch_backend("Agg")

# ---------------------------------------------------------
# Run script
# ---------------------------------------------------------

random.seed(params.seed)  # setting the seed for randomness

def fragment(sequence, version, length, location, min_length):
    """
    Fragment a given sequence.

    Parameters:
    - sequence (str): the sequence to fragment.
    - version (str): fragmentation method, either "absolute" or "percent".
    - length (int): length for fragmentation. If version is "percent", this is treated as a percentage.
    - location (str): where the fragmentation happens - "start", "end", or "random".
    - min_length (int): minimum length a fragment should be to be retained.

    Returns:
    - list of fragments.
    """
    frags = []
    if version == "absolute":
        for i in range(0, len(sequence) - length + 1):
            frags.append(sequence[i : i + length])

    elif version == "percent":
        actual_length = int(len(sequence) * (length / 100))
        for i in range(0, len(sequence) - actual_length + 1):
            frags.append(sequence[i : i + actual_length])

            # Filter fragments based on min_length
    frags = [frag for frag in frags if len(frag) >= min_length]

    # Retention logic based on the location parameter
    if location == "start":
        if frags:
            frags = [frags[0]]
    elif location == "end":
        if frags:
            frags = [frags[-1]]
    elif location == "random":
        if frags:
            chosenIndex = random.randint(0, len(frags) - 1)
            frags = [frags[chosenIndex]]

    return frags

with open(input.fasta, "r") as f:
    fastaSequences = SeqIO.parse(f, "fasta")

    with open(output.fasta_out, "w") as the_file:
        for fasta in fastaSequences:
            title_line, sequence = fasta.description, str(fasta.seq)

            fragments = fragment(
                sequence,
                params.version,
                params.frag_length,
                params.location,
                params.min_length,
            )

            for i, frag in enumerate(fragments):
                the_file.write(
                    ">" + title_line + " Fragment=" + str(i) + "\n"
                )
                the_file.write(frag + "\n")