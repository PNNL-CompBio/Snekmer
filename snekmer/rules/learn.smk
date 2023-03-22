
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
base_file = glob(join(input_dir, "base", "*.csv"))
# map extensions to basename (basename.ext.gz -> {basename: ext})
UZ_MAP = {
    skm.utils.split_file_ext(f)[0]: skm.utils.split_file_ext(f)[1] for f in zipped
}
FA_MAP = {
    skm.utils.split_file_ext(f)[0]: skm.utils.split_file_ext(f)[1] for f in unzipped
}

# get unzipped filenames
UZS = [f"{f}.{ext}" for f, ext in UZ_MAP.items()]
# isolate basenames for all files
FAS = list(FA_MAP.keys())
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
# check method
methods = ["score","score_enhanced","cosine"]
if config["learnapp"]['type'] not in methods:
    sys.exit("Please select a scoring 'type' in the config file under 'learnapp'. Options are 'score', 'score_enhanced', and 'cosine'.")
else:
    method = config["learnapp"]['type']
# check other learnapp config values that are T/F
options =[(config["learnapp"]['save_summary']),(config["learnapp"]['save_apply_associations']), (config["learnapp"]['save_results'])]
if all((option == True or option == False) for option in options) == False:
    sys.exit("Incorrect Value Selected. Please check a 'save_summary','save_apply_associations', or 'save_results' in the config file under 'learnapp'. Options are 'True' or 'False'.")
#check summary_topN
if config["learnapp"]['save_summary'] == True:
    if type(config["learnapp"]['summary_topN']) != int:
        sys.exit("Incorrect Value Selected for summary_topN. Value must be an integer greater or equal to zero.")
    if config["learnapp"]['summary_topN'] < 0:
        sys.exit("Incorrect Value Selected for summary_topN. Value must be an integer greater or equal to zero.")

# define output files to be created by snekmer
rule all:
    input:
        # expand(join(input_dir, "{uz}"), uz=UZS),  # require unzipping
        expand(join("output", "learn", "kmer-counts-{nb}.csv"), nb=FAS),
        join("output", "learn", "kmer-counts-total.csv"),
        join("output", "learn", "kmer-associations.csv"),

# if any files are gzip zipped, unzip them
use rule unzip from process with:
    output:
        unzipped=join(input_dir, "{uz}"),
        zipped=join(input_dir, "zipped", "{uz}.gz"),

# build kmer count vectors for each basis set
use rule vectorize from kmerize with:
    input:
        fasta=lambda wildcards: join(input_dir, f"{wildcards.nb}.{FA_MAP[wildcards.nb]}"),
    output:
        data=join("output", "vector", "{nb}.npz"),
        kmerobj=join("output", "kmerize", "{nb}.kmers"),
    log:
        join("output", "kmerize", "log", "{nb}.log"),
# WORKFLOW to learn kmer associations
# collect all seq files and generate mega-cluster
rule learn:
    input:
        # data=expand(join("output", "vector", "{fas}.npz"),fas = FAS),
        data="output/vector/{nb}.npz",
        #data=rules.vectorize.output.data,
        annotation=expand("{an}", an=annot_files),
        # base_counts=expand("{bf}", bf=base_file),
    output:
        table="output/learn/kmer-counts-{nb}.csv"
        # totals=join("output", "learn", "kmer-counts-total.csv"),
    log:
        join(out_dir, "learn", "log", "learn-{nb}.log"),
    run:
        # log script start time
        start_time = datetime.now()
        with open(log[0], "a") as f:
            f.write(f"start time:\t{start_time}\n")
        #Note - if the same protein (determined by seqID) is read multiple times in the same FASTA, only the last version is kept. (earlier ones are overwritten)

        annots = list()
        for f in input.annotation:
            annots.append(pd.read_table(f))
        annotations = pd.concat(annots)
        master_seq_annot_dict = {}
        master_seq_set= annots[0]['id'].tolist()
        # master_annot_set = annots[0]['InterPro'].tolist()
        master_annot_set = annots[0]['TIGRFAMs'].tolist()

        for i,seqid in enumerate(master_seq_set):
            master_seq_annot_dict[seqid] = master_annot_set[i]
        master_seq_set = set(master_seq_set)
        master_annot_set = set(master_annot_set)
        f = input.data
        # for file_num,f in enumerate(input.data):
        df, kmerlist = skm.io.load_npz(f)
        seqids = df["sequence_id"]
        # print(df)
        kmer_totals = []
        for item in kmerlist:
            kmer_totals.append(0)
        k_len = len(kmerlist[0])
        #Generate kmer counts
        seq_kmer_dict = {}
        counter = 0
        # num_of_kmers = []
        # print(seqids)
        for i,seq in enumerate(seqids):
            v = df["sequence"][i]
            kmer_counts = dict()
            items = []
            for item in range(0,(len((v)) - k_len +1)):
                items.append(v[item:(item+k_len)])
            #define in config
            if method in ["score","score_enhanced"]:
                #presence/absence
                items = set(items)
            elif method in ["cosine"]:
                #counts
                pass
            for j in items:
                kmer_counts[j] = kmer_counts.get(j, 0) + 1  
            store = []
            #convert kmer:count dict to same order as kmerlist, fill in 0s
            for i,item in enumerate(kmerlist):
                if item in kmer_counts:
                    store.append(kmer_counts[item])
                    kmer_totals[i] += kmer_counts[item]
                else:
                    store.append(0)
            #store values in kmer_seq_dict
            seq_kmer_dict[seq]= store

        #Filter out annotations not in master list
        #keep kmer counts
        #keep Sequence counts
        remove_these = []
        annotation_count_dict = {}
        total_seqs = len(seq_kmer_dict)
        # print(seq_kmer_dict)
        for i,seqid in enumerate(list(seq_kmer_dict)):
            
            x =re.findall(r'\|(.*?)\|', seqid)[0]
            if x not in master_seq_set:
                del seq_kmer_dict[seqid]
                remove_these.append(i)
            else:
                if master_seq_annot_dict[x] not in seq_kmer_dict:
                    seq_kmer_dict[master_seq_annot_dict[x]] = seq_kmer_dict.pop(seqid)
                else:
                    zipped_lists = zip(seq_kmer_dict.pop(seqid), seq_kmer_dict[master_seq_annot_dict[x]])
                    seq_kmer_dict[master_seq_annot_dict[x]] = [x + y for (x, y) in zipped_lists]
                    remove_these.append(i)
                if master_seq_annot_dict[x] not in annotation_count_dict:
                    annotation_count_dict[master_seq_annot_dict[x]] = 1
                else: 
                    annotation_count_dict[master_seq_annot_dict[x]] += 1
        #construct final output
        counts_and_sums_df = pd.DataFrame(seq_kmer_dict.values())        
        counts_and_sums_df.insert(0,"Annotations",annotation_count_dict.values(),True)
        counts_and_sums_df.insert(1,"Kmer Count",(counts_and_sums_df[list(counts_and_sums_df.columns[1:])].sum(axis=1).to_list()),True)
        kmer_totals[0:0] = [0,total_seqs]
        colnames = ["Sequence count"] + ["Kmer Count"] + list(kmerlist)
        # print(counts_and_sums_df)
        counts_and_sums_df = pd.DataFrame(np.insert(counts_and_sums_df.values, 0, values=(kmer_totals), axis=0))
        counts_and_sums_df.columns = colnames
        new_index = ["Totals"] + list(annotation_count_dict.keys())
        counts_and_sums_df.index = new_index
    
        out_name = "output/learn/kmer-counts-" + str(f)[14:-4] + ".csv"
        counts_and_sums_df_out = pa.Table.from_pandas(counts_and_sums_df,preserve_index=True)
        csv.write_csv(counts_and_sums_df_out, out_name)
        skm.utils.log_runtime(log[0], start_time)

rule merge:
    input:
        table=expand(join("output", "learn", "kmer-counts-{nb}.csv"), nb=FAS),
        base_counts=expand("{bf}", bf=base_file),
    output:
        totals=join("output", "learn", "kmer-counts-total.csv"),
    log:
        join(out_dir, "learn", "log", "merge.log"),
    run:
        # print(input.table)
        for file_num,f in enumerate(input.table):
            print("database #: ",file_num,"\n")
            counts_and_sums_df = pd.read_csv(str(f), index_col="__index_level_0__", header=0, engine="pyarrow")
            print(counts_and_sums_df)
            if file_num == 0:
                running_merge = counts_and_sums_df
            # elif file_num == 1:
            #     running_merge = (pd.concat([running_merge,counts_and_sums_df]).reset_index().groupby('__index_level_0__', sort=False).sum(min_count=1)).fillna(0)
            elif file_num >= 1:
                running_merge = (pd.concat([running_merge,counts_and_sums_df]).reset_index().groupby('__index_level_0__', sort=False).sum(min_count=1)).fillna(0)

        #Merge Step - joined with learn because it was much faster.
        base_check = False
        print("\nChecking for base file to merge with.\n")
        if "csv" in str(input.base_counts):
            print("CSV detected. Matching annotations, kmers, and totals will be summed. New annotations and kmers will be added.\n")
            base_check = True
        elif input.base_counts == "": 
            print("No base directory detected\n")
        elif str(input.base_counts) == "input/base": 
            print("Empty base directory detected\n")
        else:
            print("No file type detected. Please use a .csv file in input/base directory.\n")
        #check that kmer lengths and alphabets match base file
        if base_check == True:
            #read_csv is slow, likely replaceable
            base_df = pd.read_csv(str(input.base_counts), index_col=-1, header=0, engine="pyarrow")
            print("\nBase Database: \n")
            print(base_df)
            # Here is an assumption which is likely true but in cases with extremely high kmer values / large alphabets, the odds decrease.
            # I assume that all kmer alphabet values will be found within the first 4000 items in kmerlist.
            # This is to check that these two data structures use the same two alphabets. This is done for speed purposes, but may not be nesessary
            check_1 = 4000
            check_2 = 4000
            if len(running_merge.columns.values) < 4000:
                check_1 = len(running_merge.columns.values)
            if len(running_merge.columns.values) < 4000:
                check_1 = len(running_merge.columns.values)
            alphabet_initial = set(itertools.chain(*[list(x) for x in running_merge.columns.values[3:check_1]]))
            alphabet_base = set(itertools.chain(*[list(x) for x in base_df.columns.values[3:check_1]]))
            print(alphabet_initial)
            print(alphabet_base)
            if alphabet_base == alphabet_initial:
                base_check = True
            else: 
                base_check = False
                print("Different Alphabets Detected. Base File not merged.")
        
        if base_check == True:
            print(len(str(running_merge.columns.values[1])))
            if len(str(running_merge.columns.values[1])) == len(str(base_df.columns.values[1])):
                base_check = True
            else: 
                base_check = False
                print("Different kmer lengths detected. Base File not merged.")
        if base_check == True:
            print("\nMerged Database \n")
            xy = (pd.concat([base_df, running_merge]).reset_index().groupby('__index_level_0__', sort=False).sum(min_count=1)).fillna(0)
            xy_out = pa.Table.from_pandas(xy,preserve_index=True)
            csv.write_csv(xy_out, "output/learn/kmer-counts-total.csv")
            print(xy)
        else:
            print("\Database Merged. Not merged with base file. \n")
            running_merge_out = pa.Table.from_pandas(running_merge,preserve_index=True)
            csv.write_csv(running_merge_out, "output/learn/kmer-counts-total.csv")




rule generate_probabilities:
    input:
        totals=join("output", "learn", "kmer-counts-total.csv"),
    output:
        join("output", "learn", "kmer-associations.csv"),
    log:
        join(out_dir, "learn", "log", "learn_generate.log"),
    run:
        start_time = datetime.now()
        with open(log[0], "a") as f:
            f.write(f"start time:\t{start_time}\n")
        #read_csv is slow, maybe replaceable
        totals = pd.read_csv(str(input.totals), index_col="__index_level_0__", header=0, engine="pyarrow")
        totals
        #this line is also slow on very large datasets
        prob_2 = totals.iloc[1:, 2:].div(totals["Kmer Count"].tolist()[1:], axis = "rows")
        print(prob_2)
        prob_2 = pa.Table.from_pandas(prob_2,preserve_index=True)
        csv.write_csv(prob_2, 'output/learn/kmer-associations.csv')
        
        skm.utils.log_runtime(log[0], start_time)
