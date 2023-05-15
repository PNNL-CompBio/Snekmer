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
# Seq-Annotation-Scores
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
    

rule all:
    input:
        expand(join(input_dir, "{uz}"), uz=UZS),
        # ["output/apply/Seq-Annotation-Scores-{fa}.csv".format(fa=FAS) for dataset in FAS],
        expand(join("output","apply","Seq-Annotation-Scores-{nb}.csv"),nb=FAS),
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
        "output/apply/Seq-Annotation-Scores-{nb}.csv",
        "output/apply/kmer-summary-{nb}.csv"
    log:
        join(out_dir, "apply", "log", "{nb}.log"),

    run:
        #Note - if the same protein (determined by seqID) is read multiple times in the same FASTA, only the last version is kept. (earlier ones are overwritten)
        #read in compare association.
        #generate new association
        #check they match
        #compare them
        #check that kmer lengths and alphabets match compare file
        
            #read_csv is very slow, likely replaceable
        start_time = datetime.now()
        with open(log[0], "a") as f:
            f.write(f"start time:\t{start_time}\n")
        ####        For Validation      ########
        # annots = list()
        # for f in input.annotation:
        #     annots.append(pd.read_table(f))
        # annotations = pd.concat(annots)
        # master_seq_annot_dict = {}
        # master_seq_set= annots[0]['id'].tolist()
        # master_annot_set = annots[0]['TIGRFAMs'].tolist()
        # for i,seqid in enumerate(master_seq_set):
        #     master_seq_annot_dict[seqid] = master_annot_set[i]
        # master_seq_set = set(master_seq_set)
        # master_annot_set = set(master_annot_set)
        compare_loaded = pd.read_csv(str(input.compare_associations), index_col="__index_level_0__", header=0, engine="c")
        ########################################

        f = input.data
        print("Started:",f)
        df, kmerlist = skm.io.load_npz(f)
        seqids = df["sequence_id"]
        kmer_totals = []
        for item in kmerlist:
            kmer_totals.append(0)
        k_len = len(kmerlist[0])
        #Generate kmer counts
        seq_kmer_dict = {}
        counter = 0
        for i,seq in enumerate(seqids):
            v = df["sequence"][i]
            kmer_counts = dict()
            items = []
            #kmerize seq
            for item in range(0,(len((v)) - k_len +1)):
                items.append(v[item:(item+k_len)])
            #create kmer:count dict
            #define in config - Set Presence/Absence vs Counts
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
            
        annotation_count_dict = {}
        total_seqs = len(seq_kmer_dict)
        #construct final output
        counts_and_sums_df = pd.DataFrame(seq_kmer_dict.values())        
        counts_and_sums_df.insert(0,"Annotations",1,True)
        kmer_totals.insert(0,total_seqs)
        colnames = ["Sequence count"] + list(kmerlist)
        counts_and_sums_df = pd.DataFrame(np.insert(counts_and_sums_df.values, 0, values=kmer_totals, axis=0))
        counts_and_sums_df.columns = colnames
        new_index = ["Totals"] + list(seq_kmer_dict.keys())
        counts_and_sums_df.index = new_index
        
        #Also this particular line is slow. Can this division be done faster with linear alg...?
        new_associations = counts_and_sums_df.iloc[1:, 1:].div(counts_and_sums_df["Sequence count"].tolist()[1:], axis = "rows")
        
        #optionally write new assocations to file
        if config["learnapp"]['save_apply_associations'] == True:
            out_name = "output/apply/kmer-associations-" + str(f)[14:-4] + ".csv"
            new_associations_write = pa.Table.from_pandas(new_associations)
            csv.write_csv(new_associations_write, out_name)
        compare_df = compare_loaded.copy(deep=True)

        if method in ["score_enhanced", "score"]:
            if len(str(new_associations.columns.values[10])) == len(str(compare_df.columns.values[10])):
                compare_check = True
                # print("Same kmer lengths detected. Proceeding")
            else: 
                compare_check = False
                sys.exit("Please confirm kmer lengths match in config files of snekmer 'learn' and 'apply'.")
            if compare_check == True:
                check_1 = len(new_associations.columns.values)
                check_2 = len(compare_df.columns.values)
                alphabet_initial = set(itertools.chain(*[list(x) for x in new_associations.columns.values[10:check_1]]))
                alphabet_compare = set(itertools.chain(*[list(x) for x in compare_df.columns.values[10:check_1]]))
                if alphabet_compare == alphabet_initial:
                    compare_check = True
                else: 
                    compare_check = False

            new_cols = set(new_associations.columns)
            compare_cols = set(compare_df.columns)
            add_to_compare = []
            add_to_new = []
            for val in new_cols:
                if val not in compare_cols:
                    add_to_compare.append(val)
            for val in compare_cols:
                if val not in new_cols:
                    add_to_new.append(val)
            d = dict.fromkeys(add_to_compare, 0)
            temp_df = pd.DataFrame(d, index=compare_df.index)
            compare_df = pd.concat([compare_df, temp_df], axis=1)
            d = dict.fromkeys(add_to_new, 0)
            temp_df = pd.DataFrame(d, index=new_associations.index)
            new_associations = pd.concat([new_associations, temp_df], axis=1)

        #Cosine
        if method =="cosine":
            if len(str(counts_and_sums_df.columns.values[10])) == len(str(compare_df.columns.values[10])):
                compare_check = True
            else: 
                compare_check = False
            if compare_check == True:
                check_1 = len(new_associations.columns.values)
                check_2 = len(compare_df.columns.values)
                alphabet_initial = set(itertools.chain(*[list(x) for x in counts_and_sums_df.columns.values[10:check_1]]))
                alphabet_compare = set(itertools.chain(*[list(x) for x in compare_df.columns.values[10:check_1]]))
                if alphabet_compare == alphabet_initial:
                    compare_check = True
                    print("Same alphabets detected. Proceeding")
                else: 
                    compare_check = False
                    # print("Different Alphabets Detected. Compare File not merged.")
            new_cols = set(counts_and_sums_df.columns)
            compare_cols = set(compare_df.columns)
            add_to_compare = []
            add_to_new = []
            for val in new_cols:
                if val not in compare_cols:
                    add_to_compare.append(val)
            for val in compare_cols:
                if val not in new_cols:
                    add_to_new.append(val)
            d = dict.fromkeys(add_to_compare, 0)
            temp_df = pd.DataFrame(d, index=compare_df.index)
            compare_df = pd.concat([compare_df, temp_df], axis=1)
            d = dict.fromkeys(add_to_new, 0)
            temp_df = pd.DataFrame(d, index=counts_and_sums_df.index)
            counts_and_sums_df = pd.concat([counts_and_sums_df, temp_df], axis=1)

        
        #Cosine similarity Score
        #A note to user. 
        #Cosine Similarity score requires COUNTS data
        #Score methods require ASSOCIATION data
        if method == "cosine":
            #remove undesired rows and columns
            counts_and_sums_df.drop(columns=counts_and_sums_df.columns[-1:], axis=1,  inplace=True)
            counts_and_sums_df.drop(columns=counts_and_sums_df.columns[:1], axis=1,  inplace=True)
            compare_df.drop(columns=compare_df.columns[:2], axis=1,  inplace=True)
            counts_and_sums_df.drop(index="Totals", axis=0, inplace=True)
            compare_df.drop(index="Totals", axis=0, inplace=True)
            d = sklearn.metrics.pairwise.cosine_similarity(compare_df,counts_and_sums_df).T
            final_matrix_with_scores = pd.DataFrame(d, columns=compare_df.index, index=counts_and_sums_df.index)
            counts_and_sums_df = []
            compare_df =[]
        if method == "score_enhanced":
            compare_df = compare_df.replace(0,np.nan)
            compare_df = compare_df.where(compare_df.notna(), -compare_df.mean(axis=1), axis=0)
        if method == "score" or method == "score_enhanced":
            d=(new_associations.values@compare_df.values.T)
            final_matrix_with_scores = pd.DataFrame(d, columns=compare_df.index, index=new_associations.index)
        if method == "experimental":
            #create your own here!
            #manhattan distance - least accurate 
            # counts_and_sums_df = counts_and_sums_df.drop('Sequence count', axis=1)
            # counts_and_sums_df = counts_and_sums_df.drop('Totals', axis=0)
            # d = sklearn.metrics.pairwise.manhattan_distances(compare_df,counts_and_sums_df).T
            #note, you will need to choose  Presence/Absence or Counts above (and in Learn).
            pass

        #write output
        if config["learnapp"]['save_results'] == True: 
            out_name = "output/apply/Seq-Annotation-Scores-" + str(f)[14:-4] + ".csv"
            final_matrix_with_scores_write = pa.Table.from_pandas(final_matrix_with_scores)
            csv.write_csv(final_matrix_with_scores_write, out_name)

        ########
        #create output
        # Protein ID, Prediction, Score, delta, Confidence
        global_conf_scoring = pd.read_csv(str(input.confidence_associations))
        global_conf_scoring.index= global_conf_scoring[global_conf_scoring.columns[0]]
        global_conf_scoring = global_conf_scoring.iloc[: , 1:]
        global_conf_scoring = global_conf_scoring[global_conf_scoring.columns[0]].squeeze()

        print("global_conf_scoring: \n",global_conf_scoring)

        w =[]
        order = np.argsort(-final_matrix_with_scores.values, axis=1)[:, :2]
        print("order: ", order)
        for i,item in enumerate(order):
            w.append((final_matrix_with_scores[final_matrix_with_scores.columns[[item]]][i:i+1]).values.tolist()[0])
            
        diffs = []
        top1 = []
        for item in w:
            diffs.append(item[0] - item[1])
            top1.append(item[0])
    
        vals = pd.DataFrame({'delta':diffs})

        names = pd.DataFrame(final_matrix_with_scores.columns[order][:, :1])
        score = pd.DataFrame(top1)
        score.columns = ["Score"]
        names.columns = ["Prediction"]
        names = names.astype(str)

        vals = vals.round(decimals=2)
        #map value in vals to score
        vals['Confidence'] = vals["delta"].map(global_conf_scoring)

        
        summary = pd.concat([names,score,vals],axis=1)
        summary.index=final_matrix_with_scores.index

        out_name_2 = "output/apply/kmer-summary-" + str(f)[14:-4] + ".csv"

        summary.reset_index(inplace=True)
        summary_write = pa.Table.from_pandas(summary)
        csv.write_csv(summary_write, out_name_2)

        skm.utils.log_runtime(log[0], start_time)