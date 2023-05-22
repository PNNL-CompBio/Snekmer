# force snakemake v6.0+ (required for modules)
from snakemake.utils import min_version
min_version("6.0")
# load snakemake modules
module process:
    snakefile:
        "process.smk"
    config:
        config
module kmerize:
    snakefile:
        "kmerize.smk"
    config:
        config
# built-in imports
import gzip
import json
import pickle
import struct
from datetime import datetime
from glob import glob
from itertools import product, repeat
from multiprocessing import Pool
from os import makedirs
from os.path import basename, dirname, exists, join, splitext, split
import numpy as np
import pandas as pd
import seaborn as sns
import snekmer as skm
from Bio import SeqIO
from pathlib import Path
import copy
from scipy.stats import rankdata
import csv
import sys
import time
import pyarrow as pa
import pyarrow.csv as csv
import itertools
from itertools import islice
from scipy import spatial
import sys
import sklearn
#Note:
#Pyarrow installed via "conda install -c conda-forge pyarrow"
# collect all fasta-like files, unzipped filenames, and basenames
input_dir = "input" if (("input_dir" not in config) or (str(config["input_dir"]) == "None")) else config["input_dir"]
input_files = glob(join(input_dir, "*"))
# base_file = glob(join(input_dir,"base" "*"))
zipped = [fa for fa in input_files if fa.endswith(".gz")]
unzipped = [
    fa.rstrip(".gz")
    for fa, ext in product(input_files, config["input_file_exts"])
    if fa.rstrip(".gz").endswith(f".{ext}")
]
annot_files = glob(join("annotations", "*.ann"))
compare_file = glob(join("counts", "*.csv"))
confidence_file = glob(join("confidence", "*.csv"))
# map extensions to basename (basename.ext.gz -> {basename: ext})
UZ_MAP = {
    skm.utils.split_file_ext(f)[0]: skm.utils.split_file_ext(f)[1] for f in zipped
}
FA_MAP = {
    skm.utils.split_file_ext(f)[0]: skm.utils.split_file_ext(f)[1] for f in unzipped
}
# seq-annotation-scores
# get unzipped filenames
UZS = [f"{f}.{ext}" for f, ext in UZ_MAP.items()]
# isolate basenames for all files
FAS = list(FA_MAP.keys())
FAV = list(FA_MAP.values())
# parse any background files
bg_files = glob(join(input_dir, "background", "*"))
if len(bg_files) > 0:
    bg_files = [skm.utils.split_file_ext(basename(f))[0] for f in bg_files]
NON_BGS, BGS = [f for f in FAS if f not in bg_files], bg_files
# terminate with error if invalid alphabet specified
skm.alphabet.check_valid(config["alphabet"])
# define output directory (helpful for multiple runs)
out_dir = skm.io.define_output_dir(
    config["alphabet"], config["k"], nested=config["nested_output"]
)

wildcard_constraints:
    dataset=FAS,
    FAS=FAS


options =[(config["learnapp"]['save_apply_associations'])]
if all((option == True or option == False) for option in options) == False:
    sys.exit("Incorrect Value Selected. Please check if 'save_apply_associations' in in the config file under 'learnapp'. Options are 'True' or 'False'.")


rule all:
    input:
        expand(join(input_dir, "{uz}"), uz=UZS),
        # ["output/apply/seq-annotation-scores-{fa}.csv".format(fa=FAS) for dataset in FAS],
        expand(join("output","apply","seq-annotation-scores-{nb}.csv"),nb=FAS),
        expand(join("output","apply","kmer-summary-{nb}.csv"),nb=FAS),

use rule vectorize from kmerize with:
    input:
        fasta=lambda wildcards: join("input", f"{wildcards.nb}.{FA_MAP[wildcards.nb]}"),
    output:
        data=join("output", "vector", "{nb}.npz"),
        kmerobj=join("output", "kmerize", "{nb}.kmers"),
    log:
        join("output", "kmerize", "log", "{nb}.log"),

rule apply:
    input:
        data="output/vector/{nb}.npz",
        annotation=expand("{an}", an=annot_files),
        compare_associations=expand("{comp}", comp=compare_file),
        confidence_associations=expand("{conf}", conf=confidence_file),
    output:
        "output/apply/seq-annotation-scores-{nb}.csv",
        "output/apply/kmer-summary-{nb}.csv"
    log:
        join(out_dir, "apply", "log", "{nb}.log"),

    run:
        start_time = datetime.now()
        with open(log[0], "a") as f:
            f.write(f"start time:\t{start_time}\n")

        ##### Generate Inputs
        kmer_count_totals = pd.read_csv(str(input.compare_associations), index_col="__index_level_0__", header=0, engine="c")
        df, kmerlist = skm.io.load_npz(input.data)
        seqids = df["sequence_id"]
        kmer_totals = []
        for item in kmerlist:
            kmer_totals.append(0)
        k_len = len(kmerlist[0])

        ##### Generate Kmer Counts
        seq_kmer_dict = {}
        counter = 0
        for i,seq in enumerate(seqids):
            v = df["sequence"][i]
            kmer_counts = dict()
            items = []
            for item in range(0,(len((v)) - k_len +1)):
                items.append(v[item:(item+k_len)])
            for j in items:
                kmer_counts[j] = kmer_counts.get(j, 0) + 1  
            store = []
            for i,item in enumerate(kmerlist):
                if item in kmer_counts:
                    store.append(kmer_counts[item])
                    kmer_totals[i] += kmer_counts[item]
                else:
                    store.append(0)
            seq_kmer_dict[seq]= store
            
        
        ######  Construct Kmer Counts Dataframe
        total_seqs = len(seq_kmer_dict)
        kmer_counts = pd.DataFrame(seq_kmer_dict.values())        
        kmer_counts.insert(0,"Annotations",1,True)
        kmer_totals.insert(0,total_seqs)
        kmer_counts = pd.DataFrame(np.insert(kmer_counts.values, 0, values=kmer_totals, axis=0))
        kmer_counts.columns = ["Sequence count"] + list(kmerlist)
        kmer_counts.index = ["Totals"] + list(seq_kmer_dict.keys())
        
        new_associations = kmer_counts.iloc[1:, 1:].div(kmer_counts["Sequence count"].tolist()[1:], axis = "rows")
        
        ##### Make Kmer Counts Dataframe match Kmer Counts Totals Format
        if len(str(kmer_counts.columns.values[10])) == len(str(kmer_count_totals.columns.values[10])):
            compare_check = True
        else: 
            compare_check = False
        if compare_check == True:
            check_1 = len(new_associations.columns.values)
            check_2 = len(kmer_count_totals.columns.values)
            alphabet_initial = set(itertools.chain(*[list(x) for x in kmer_counts.columns.values[10:check_1]]))
            alphabet_compare = set(itertools.chain(*[list(x) for x in kmer_count_totals.columns.values[10:check_1]]))
            if alphabet_compare == alphabet_initial:
                compare_check = True
        if compare_check == False:
            print("Compare Check Failed. ")
            sys.exit()

        new_cols = set(kmer_counts.columns)
        compare_cols = set(kmer_count_totals.columns)
        add_to_compare = []
        add_to_new = []
        for val in new_cols:
            if val not in compare_cols:
                add_to_compare.append(val)
        for val in compare_cols:
            if val not in new_cols:
                add_to_new.append(val)

        kmer_count_totals = pd.concat([kmer_count_totals, pd.DataFrame(dict.fromkeys(add_to_compare, 0), index=kmer_count_totals.index)], axis=1)
        kmer_count_totals.drop(columns=kmer_count_totals.columns[:2], index="Totals", axis=0, inplace=True)
        kmer_counts = pd.concat([kmer_counts, pd.DataFrame(dict.fromkeys(add_to_new, 0), index=kmer_counts.index)], axis=1)
        kmer_counts.drop(columns=kmer_counts.columns[-1:].union(kmer_counts.columns[:1]), index="Totals", axis=0, inplace=True)


        #### Perform Cosine Similarity between Kmer Counts Totals and Counts and Sums DF
        cosine_df = sklearn.metrics.pairwise.cosine_similarity(kmer_count_totals,kmer_counts).T
        kmer_count_totals = pd.DataFrame(cosine_df, columns=kmer_count_totals.index, index=kmer_counts.index)

        ##### Write Optional Output
        if config["learnapp"]['save_results'] == True: 
            out_name = "output/apply/seq-annotation-scores-" + str(input.data)[14:-4] + ".csv"
            kmer_count_totals_write = pa.Table.from_pandas(kmer_count_totals)
            csv.write_csv(kmer_count_totals_write, out_name)

        
        ##### Create True Output
        # Protein ID, Prediction, Score, delta, Confidence
        global_confidence_scores = pd.read_csv(str(input.confidence_associations))
        global_confidence_scores.index= global_confidence_scores[global_confidence_scores.columns[0]]
        global_confidence_scores = global_confidence_scores.iloc[: , 1:]
        global_confidence_scores = global_confidence_scores[global_confidence_scores.columns[0]].squeeze()

        score_rank =[]
        sorted_vals = np.argsort(-kmer_count_totals.values, axis=1)[:, :2]
        for i,item in enumerate(sorted_vals):
            score_rank.append((kmer_count_totals[kmer_count_totals.columns[[item]]][i:i+1]).values.tolist()[0])
            
        delta = []
        Top_Score = []
        for score in score_rank:
            delta.append(score[0] - score[1])
            Top_Score.append(score[0])
    
        vals = pd.DataFrame({'delta':delta})
        predictions = pd.DataFrame(kmer_count_totals.columns[sorted_vals][:, :1])
        score = pd.DataFrame(Top_Score)
        score.columns = ["Score"]
        predictions.columns = ["Prediction"]
        predictions = predictions.astype(str)
        vals = vals.round(decimals=2)
        vals['Confidence'] = vals["delta"].map(global_confidence_scores)

        results = pd.concat([predictions,score,vals],axis=1)
        results.index=kmer_count_totals.index

        #### Write Results 
        out_name_2 = "output/apply/kmer-summary-" + str(input.data)[14:-4] + ".csv"
        results.reset_index(inplace=True)
        results_write = pa.Table.from_pandas(results)
        csv.write_csv(results_write, out_name_2)

        skm.utils.log_runtime(log[0], start_time)