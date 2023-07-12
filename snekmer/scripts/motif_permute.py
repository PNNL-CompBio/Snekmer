"""
Created on Tue Apr 25 15:14:39 2023

author: @tnitka
"""

# ---------------------------------------------------------
# Imports
# ---------------------------------------------------------
import pickle

import snekmer as skm
import pandas as pd
import gzip

# ---------------------------------------------------------
# Files and parameters
# ---------------------------------------------------------

config = snakemake.config

# ---------------------------------------------------------
# Run script
# ---------------------------------------------------------
    
# load input data
with open(snakemake.input.matrix, "rb") as f:
    data = pickle.load(f)
    
# with open(snakemake.input.kmers, "rb") as f:
#     kmers = f.readlines()
    
with gzip.open(snakemake.input.weights, "rb") as f:
    weights = pd.read_csv(f)
    
# prevent kmer NA being read as np.nan
if config["k"] == 2:
    weights["kmer"] = weights["kmer"].fillna("NA")

kmers = weights['kmer'].values    
scores = weights['sample'].values
family = skm.utils.get_family(
    skm.utils.split_file_ext(snakemake.input.weights)[0],
    regex=config["input_file_regex"],
)
scorer = skm.score.KmerScorer()

# set number of permutations to test
n_iter = (
    config["motif"]["n"]  
    )


# get kmers for this particular set of sequences
with open(snakemake.input.kmerobj, "rb") as f:
    kmerobj = pickle.load(f)
    
# set category label name (e.g. "family")
label = config["score"]["lname"] if str(config["score"]["lname"]) != "None" else "label"

# binary T/F for classification into family
family = skm.utils.get_family(snakemake.wildcards.nb)
binary_labels = [True if value == family else False for value in data[label]]

# get alphabet name
if config["alphabet"] in skm.alphabet.ALPHABET_ORDER.keys():
    alphabet_name = skm.alphabet.ALPHABET_ORDER[config["alphabet"]].capitalize()
else:
    alphabet_name = str(config["alphabet"]).capitalize()

  
# permute input data

score_matrix = pd.DataFrame({'kmer': kmers})
score_array = pd.DataFrame.to_numpy(score_matrix)
motif = skm.motif.SnekmerMotif()
perm_data = motif.permute(
    data, 
    skm.utils.get_family(snakemake.wildcards.nb, regex=config["input_file_regex"]),
    label_col=label)

    
# save output
perm_data.to_csv(snakemake.output.data, index=False, compression="gzip")

# record script endtime
#skm.utils.log_runtime(snakemake.log[0], start_time)