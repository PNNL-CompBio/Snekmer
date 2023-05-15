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
from scipy.interpolate import interp1d

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
compare_file = glob(join(input_dir, "compare", "*.csv"))
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
    
# define output files to be created by snekmer
# rule all:
#     input:
#         # expand(join(input_dir, "{uz}"), uz=UZS),  # require unzipping
#         "output/apply/Seq-Annotation-Scores-{FAS}.csv".format(FAS=FAS),
#         # ["output/apply/Seq-Annotation-Scores-{dataset}.csv".format(dataset=dataset) for dataset in FAS],



rule all:
    input:
        expand(join(input_dir, "{uz}"), uz=UZS),
        # ["output/apply/Seq-Annotation-Scores-{fa}.csv".format(fa=FAS) for dataset in FAS],
        expand(join("output","eval_apply","Seq-Annotation-Scores-{nb}.csv"),nb=FAS),
        # expand(join("output","eval_conf","confidence_matrix-{nb}.csv"),nb=FAS),
        "output/eval_conf/confidence_matrix.csv",

use rule vectorize from kmerize with:
    input:
        fasta=lambda wildcards: join("input", f"{wildcards.nb}.{FA_MAP[wildcards.nb]}"),
    output:
        data=join("output", "vector", "{nb}.npz"),
        kmerobj=join("output", "kmerize", "{nb}.kmers"),
    log:
        join("output", "kmerize", "log", "{nb}.log"),
    # params:
    #     FAS=lambda wildcards: wildcards.FAS,
        # FAV=lambda wildcards: wildcards.FAV,

rule eval_apply:
    input:
        data="output/vector/{nb}.npz",
        annotation=expand("{an}", an=annot_files),
        compare_associations=expand("{bf}", bf=compare_file),
    output:
        "output/eval_apply/Seq-Annotation-Scores-{nb}.csv"
        # 'output/eval_apply/Complete.txt'
    log:
        join(out_dir, "eval_apply", "log", "{nb}.log"),
    # params:
    #     FAS=lambda wildcards: FAS
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
        startTime = time.time()
        # print("DATA:",input.data)
        ####        For Validation      ########
        annots = list()
        for f in input.annotation:
            annots.append(pd.read_table(f))
        annotations = pd.concat(annots)
        master_seq_annot_dict = {}
        master_seq_set= annots[0]['id'].tolist()
        master_annot_set = annots[0]['TIGRFAMs'].tolist()
        for i,seqid in enumerate(master_seq_set):
            master_seq_annot_dict[seqid] = master_annot_set[i]
        master_seq_set = set(master_seq_set)
        master_annot_set = set(master_annot_set)
        #moved up here now
        compare_loaded = pd.read_csv(str(input.compare_associations), index_col="__index_level_0__", header=0, engine="c")
        ########################################
        # for f in input.data:
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
        #######     For Validation      ######
        #this adds _known and _unknown depending on if its a known annotation. This allows us to do validations and RUO.
        annotation_count_dict = {}
        total_seqs = len(seq_kmer_dict)
        count = 0
        for seqid in list(seq_kmer_dict):
            x =re.findall(r'\|(.*?)\|', seqid)[0]
            if x not in master_seq_set:
                seq_kmer_dict[(x+"_unknown_"+str(count))] = seq_kmer_dict.pop(seqid)
            else: 
                seq_kmer_dict[(master_seq_annot_dict[x] + "_known_" + str(count))] = seq_kmer_dict.pop(seqid)
            count +=1
        #######
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
        executionTime = (time.time() - startTime)
        
        #Also this particular line is slow. Can this division be done faster with linear alg...?
        new_associations = counts_and_sums_df.iloc[1:, 1:].div(counts_and_sums_df["Sequence count"].tolist()[1:], axis = "rows")
        
        #optionally write new assocations to file
        if config["learnapp"]['save_apply_associations'] == True:
            out_name = "output/eval_apply/kmer-associations-" + str(f)[14:-4] + ".csv"
            new_associations_write = pa.Table.from_pandas(new_associations)
            csv.write_csv(new_associations_write, out_name)
        startTime = time.time()
        compare_df = compare_loaded.copy(deep=True)
        executionTime = (time.time() - startTime)
        startTime = time.time()
        if method in ["score_enhanced", "score"]:
            if len(str(new_associations.columns.values[10])) == len(str(compare_df.columns.values[10])):
                compare_check = True
                # print("Same kmer lengths detected. Proceeding")
            else: 
                compare_check = False
                # print("Different kmer lengths detected. Compare File and ", str(f)[14:-4], " don't match.")
                sys.exit("Please confirm kmer lengths match in config files of snekmer 'learn' and 'apply'.")
                # Here is an assumption which is likely true but in cases with extremely high kmer values / large alphabets, the odds decrease.
                # I assume that all kmer alphabet values will be found within the first 4000 items in kmerlist.
                # This is to check that these two data structures use the same two alphabets. 
                # it could just be set to length but would take a bit longer...
            if compare_check == True:
                check_1 = 4000
                check_2 = 4000
                if len(new_associations.columns.values) < check_1:
                    check_1 = len(new_associations.columns.values)
                if len(new_associations.columns.values) < check_2:
                    check_1 = len(new_associations.columns.values)
                alphabet_initial = set(itertools.chain(*[list(x) for x in new_associations.columns.values[10:check_1]]))
                alphabet_compare = set(itertools.chain(*[list(x) for x in compare_df.columns.values[10:check_1]]))
                if alphabet_compare == alphabet_initial:
                    compare_check = True
                else: 
                    compare_check = False
                    # print("Different Alphabets Detected. Compare File not merged.")
            executionTime = (time.time() - startTime)
            # print(' make sure files match - Execution time in seconds: ' + str(executionTime))
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
                # print("Different kmer lengths detected. Compare File and ", str(f)[14:-4], " don't match.")
                # Here is an assumption which is likely true but in cases with extremely high kmer values / large alphabets, the odds decrease.
                # I assume that all kmer alphabet values will be found within the first 4000 items in kmerlist.
                # This is to check that these two data structures use the same two alphabets. 
                # it could just be set to length but would take a bit longer...
            if compare_check == True:
                check_1 = 4000
                check_2 = 4000
                if len(counts_and_sums_df.columns.values) <= check_1:
                    check_1 = (len(counts_and_sums_df.columns.values) -5)
                if len(counts_and_sums_df.columns.values) <= check_2:
                    check_1 = (len(counts_and_sums_df.columns.values) -5)
                alphabet_initial = set(itertools.chain(*[list(x) for x in counts_and_sums_df.columns.values[10:check_1]]))
                alphabet_compare = set(itertools.chain(*[list(x) for x in compare_df.columns.values[10:check_1]]))
                if alphabet_compare == alphabet_initial:
                    compare_check = True
                    # print("Same alphabets detected. Proceeding")
                else: 
                    compare_check = False
                    # print("Different Alphabets Detected. Compare File not merged.")
            executionTime = (time.time() - startTime)
            # print(' make sure files match - Execution time in seconds: ' + str(executionTime))
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

        executionTime = (time.time() - startTime)
        startTime = time.time()
        
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

            print("long calc started for: ", f)
            #linalg
            # d=(counts_and_sums_df.values@compare_df.values.T) / np.linalg.norm(counts_and_sums_df.values,axis=1).reshape(-1,1) / np.linalg.norm(compare_df.values, axis=1)
            #einsum
            # d=np.einsum('ij,kj', counts_and_sums_df.values, compare_df.values) / np.linalg.norm(counts_and_sums_df.values,axis=1).reshape(-1,1) / np.linalg.norm(compare_df.values, axis=1)
            #sklearn - fastest for large datasets
            d = sklearn.metrics.pairwise.cosine_similarity(compare_df,counts_and_sums_df).T
            print("long calc finished for: ", f)
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
        executionTime = (time.time() - startTime)    
        # print('rest of slow time: ' + str(executionTime))
        startTime = time.time()
        #write output
        if config["learnapp"]['save_results'] == True: 
            out_name = "output/eval_apply/Seq-Annotation-Scores-" + str(f)[14:-4] + ".csv"
            final_matrix_with_scores_write = pa.Table.from_pandas(final_matrix_with_scores)
            csv.write_csv(final_matrix_with_scores_write, out_name)

        #Create summary with top N Values
        #####
        if config["learnapp"]["save_summary"] == True:
            summary_topN = config["learnapp"]['summary_topN']
            nlargest = summary_topN
            if len(final_matrix_with_scores.axes[1]) < summary_topN:
                nlargest = len(final_matrix_with_scores.axes[1])
            order = np.argsort(-final_matrix_with_scores.values, axis=1)[:, :nlargest]
            w =[]
            for i,item in enumerate(order):
                w.append((final_matrix_with_scores[final_matrix_with_scores.columns[[item]]][i:i+1]).values.tolist()[0])
            vals = pd.DataFrame(w)
            names = pd.DataFrame(final_matrix_with_scores.columns[order])
            vals = vals.astype(str)
            names = names.astype(str)
            summary = pd.concat([names[col] + "," + vals[col] for col in names],axis="columns")
            summary.columns = ['top{}'.format(i) for i in range(1, nlargest+1)]
            summary.index=final_matrix_with_scores.index
            
            # print(summary)
            out_name_2 = "output/eval_apply/kmer-summary-" + str(f)[14:-4] + ".csv"
            summary_write = pa.Table.from_pandas(summary)
            csv.write_csv(summary_write, out_name_2)
            final_matrix_with_scores =[]
            print("Completed: ",f)
            # ######
        skm.utils.log_runtime(log[0], start_time)




rule evaluate:
    input:
        # data="output/vector/{nb}.npz",
        # annotation=expand("{an}", an=annot_files),
        # compare_associations=expand("{bf}", bf=compare_file),
        # eval_apply_data="output/eval_apply/Seq-Annotation-Scores-{nb}.csv",
        # eval_apply_data= "output/eval_apply"
        eval_apply_data=expand(join("output", "eval_apply", "Seq-Annotation-Scores-{nb}.csv"),nb = FAS),
    output:
        #"output/eval_apply/Seq-Annotation-Scores-{nb}.csv"
        # 'output/eval_apply/Complete.txt'
        "output/eval_conf/confidence_matrix.csv",
    # log:
    #     join(out_dir, "eval_conf", "log", "conf.log"),
    # params:
    #     FAS=lambda wildcards: FAS
    run:
        
            #read_csv is very slow, likely replaceable
        # start_time = datetime.now()
        # with open(log[0], "a") as f:
        #     f.write(f"start time:\t{start_time}\n")
        # startTime = time.time()

        # print("aaaaa its: ", input.eval_apply_data)

        for j,f in enumerate(input.eval_apply_data):
            data = pd.read_csv(f, index_col="__index_level_0__", header=0, engine="c")
            # print(data)

            maxValueIndex = data.idxmax(axis="columns")
            result = maxValueIndex.keys()
            # print(maxValueIndex)
            # print(maxValueIndex[0])
            # print(result[0])
            # print(maxValueIndex[1])
            # print(result[1])
            

            discovered_fams = list()
            for i,item in enumerate(list(maxValueIndex)):
                discovered_fams.append(item)
            discovered_fams = set(discovered_fams)


            TF = list()
            for i,item in enumerate(list(maxValueIndex)):
                if item in result[i]:
                    TF.append("T")
                else:
                    TF.append("F")

            Known = list()
            for i,item in enumerate(list(maxValueIndex)):
                if "unknown" in result[i]:
                    Known.append("Unknown")
                else:
                    Known.append("Known")


            a = data.values
            x = a[np.arange(len(data))[:,None],np.argpartition(-a,np.arange(2),axis=1)[:,:2]]
            # print(x)
            # print(np.diff(x, axis=1))
            diff_df = pd.DataFrame(x, columns = ['Top','Second'])
            diff_df['Difference'] = -(np.diff(x, axis=1).round(decimals=2))
            diff_df['Prediction'] = list(maxValueIndex)
            diff_df['Actual'] = result
            diff_df["T/F"] = TF
            diff_df["Known/Unknown"] = Known

            known_true_diff_df = diff_df[(diff_df["Known/Unknown"] == "Known") & (diff_df["T/F"] == "T")]

            known_false_diff_df = diff_df[(diff_df["Known/Unknown"] == "Known") & (diff_df["T/F"] == "F")]

            # print(known_true_diff_df)
            # print(known_false_diff_df)
            


            possible_vals = [round(x * 0.01,2) for x in range(0, 101)]
            possible_fams = discovered_fams

            # print("possible_vals: ", possible_vals)
            # print("possible_fams: ", possible_fams)



            true_crosstab = pd.crosstab(known_true_diff_df.Prediction,known_true_diff_df.Difference)
            false_crosstab = pd.crosstab(known_false_diff_df.Prediction,known_false_diff_df.Difference)

            if j == 0:
                true_running_crosstab = true_crosstab
                false_running_crosstab = false_crosstab
            else:
                # running_crosstab = (pd.concat([running_crosstab,crosstab])).fillna(0)
                true_running_crosstab = (pd.concat([true_running_crosstab,true_crosstab]).reset_index().groupby('Prediction', sort=False).sum(min_count=1)).fillna(0)
                false_running_crosstab = (pd.concat([false_running_crosstab,false_crosstab]).reset_index().groupby('Prediction', sort=False).sum(min_count=1)).fillna(0)



            print(possible_vals)

            print(true_running_crosstab.columns)

            # true_running_crosstab.index = true_running_crosstab.index.astype('float64')
            # false_running_crosstab.index = false_running_crosstab.index.astype('float64')

            true_overlap = sorted(list(set(possible_vals) & set(true_running_crosstab.columns)))
            true_non_overlap = sorted(list(set(possible_vals) ^ set(true_running_crosstab.columns.astype(float))))
            false_overlap = sorted(list(set(possible_vals) & set(false_running_crosstab.columns)))
            false_non_overlap = sorted(list(set(possible_vals) ^ set(false_running_crosstab.columns.astype(float))))

            true_non_overlap_str = list(map(str, true_non_overlap))
            true_non_overlap_dict = dict.fromkeys(true_non_overlap_str,0)

            false_non_overlap_str = list(map(str, false_non_overlap))
            false_non_overlap_dict = dict.fromkeys(false_non_overlap_str,0)



            add_to_false = sorted(set(true_running_crosstab.index) - set(false_running_crosstab.index))

            add_to_true = sorted(set(false_running_crosstab.index) - set(true_running_crosstab.index))
            



            dex = ["TIGR0","TIGR1","TIGR2"]
            add_to_true_df = pd.DataFrame(0, index = add_to_true, columns= true_running_crosstab.columns)
            add_to_false_df = pd.DataFrame(0, index = add_to_false, columns= false_running_crosstab.columns)
            
            true_running_crosstab = pd.concat([true_running_crosstab,add_to_true_df])
            false_running_crosstab = pd.concat([false_running_crosstab,add_to_false_df])



            # print("true_overlap: ", true_overlap)
            # print("true_non_overlap: ", true_non_overlap)
            # print("false_overlap: ",false_overlap)
            # print("False_non_overlap: ", false_non_overlap)


            

            # print("true_crosstab: ", true_running_crosstab)
            # print("true_crosstab columns: ", true_running_crosstab.columns)
            # print("true_crosstab type: ", type(true_running_crosstab))

            true_running_crosstab = true_running_crosstab[true_overlap]
            # true_running_crosstab = true_running_crosstab[true_non_overlap_str] = 0
            true_running_crosstab = true_running_crosstab.assign(**true_non_overlap_dict)




            false_running_crosstab = false_running_crosstab[false_overlap]
            # false_running_crosstab = false_running_crosstab[false_non_overlap_str] = 0
            false_running_crosstab = false_running_crosstab.assign(**false_non_overlap_dict)



            #to do:

            #add columns of zeros to the true and false crosstabs.
            #these should be the missing columns (ie, the ones that dont match to the other)
            #that will allow the calculation below to proceed without resulting in NAs.


            #sort by index
            true_running_crosstab.index.names = ['Prediction']
            false_running_crosstab.index.names = ['Prediction']
            true_running_crosstab.sort_index(inplace=True) 
            false_running_crosstab.sort_index(inplace=True) 

            #sort by column
            true_running_crosstab.columns = true_running_crosstab.columns.astype(float)
            false_running_crosstab.columns = false_running_crosstab.columns.astype(float)

            true_running_crosstab = true_running_crosstab[sorted(true_running_crosstab.columns)]
            false_running_crosstab = false_running_crosstab[sorted(false_running_crosstab.columns)]


            print("true_crosstab: ", true_running_crosstab)
            print("false_crosstab: ", false_running_crosstab)

            #I think ratio could just be calculated once. un tab this.
            ratio_running_crosstab = (true_running_crosstab/(true_running_crosstab + false_running_crosstab))


            print("ratio_crosstab: ", ratio_running_crosstab)
            print("Dataframes joined: ", f, " out of ",len(input.eval_apply_data) , ".")


        true_total_dist = true_running_crosstab.sum(numeric_only=True, axis=0)
        false_total_dist = false_running_crosstab.sum(numeric_only=True, axis=0)



        ratio_total_dist = (true_running_crosstab.sum(numeric_only=True, axis=0)/(true_running_crosstab.sum(numeric_only=True, axis=0) + false_running_crosstab.sum(numeric_only=True, axis=0))).fillna(1)


        # print(true_total_dist)
        with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified also
            print("True:", true_total_dist)

        with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified also
            print("False:", false_total_dist)
        
        with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified also
            print("Ratio:", ratio_total_dist)


        print("Ratio DF: \n", ratio_running_crosstab)


        csv.write_csv(pa.Table.from_pandas(true_running_crosstab), "output/eval_conf/true_total.csv")
        csv.write_csv(pa.Table.from_pandas(false_running_crosstab), "output/eval_conf/false_total.csv")
        csv.write_csv(pa.Table.from_pandas(ratio_running_crosstab), "output/eval_conf/ratio_total.csv")


        # csv.write_csv(pa.Table.from_pandas(ratio_running_crosstab), "output/eval_conf/confidence_matrix.csv")

        def nan_helper(row):
            return np.isnan(row), lambda z: z.nonzero()[0]

        #this may mess up alignment
        #thresh can be used to set number of non-nas required, remove how
        ratio_running_crosstab = ratio_running_crosstab.dropna(axis=0,how="all")


        ratio_np_df = ratio_running_crosstab.to_numpy()
        imputed_df = np.empty((0,len(ratio_running_crosstab.columns)))

        print("imputed_df before:", imputed_df)
        print("length ratio columns", len(ratio_running_crosstab.columns))

        for row in ratio_np_df:
            nans, x= nan_helper(row)
            # print("row before", row)
            row[nans]= np.interp(x(nans), x(~nans), row[~nans])
            # print("row after", row)

            imputed_df = np.append(imputed_df, [np.array(row)], axis=0)

        imputed_df = pd.DataFrame(imputed_df,index=ratio_running_crosstab.index, columns=ratio_running_crosstab.columns)
        print("Imputed DF: \n", imputed_df)

        imputed_df = pa.Table.from_pandas(imputed_df)
        
        csv.write_csv(imputed_df, "output/eval_conf/confidence_matrix_imputed.csv")
        csv.write_csv(pa.Table.from_pandas(ratio_running_crosstab), "output/eval_conf/confidence_matrix.csv")







        # ratio_running_crosstab = ratio_running_crosstab.replace({np.nan:None})
        # x = [float(c) for c in ratio_running_crosstab.columns]
        # y = ratio_running_crosstab.values

        # f = interp1d(x, y, axis=1, kind='quadratic', fill_value='extrapolate')
        # interpolated_values = f(x)
        # interpolated_df = pd.DataFrame(interpolated_values, index=ratio_running_crosstab.index, columns=ratio_running_crosstab.columns)

        # print(interpolated_df)




            # stop here

            #possible counts table
            #https://stackoverflow.com/questions/71306034/create-count-matrix-based-on-two-columns-and-respective-range-values-pandas

            #create Two matrices that are:
            # True Positives & False Positives
            #Then join to find the ratio.
            #            .00 .01 .02 .03 .04 .05 .06 .07 .08 .09 .10  ...      (score diff)
            #  family1   
            #  family2
            #  family3
            #  .
            #  .
            #  .



            

            


            # diff_df['Correctly_Predicted'] = np.diff(x, axis=1).tolist()
            # print(diff_df)
            
























        # skm.utils.log_runtime(log[0], start_time)













# Determine if each protein annotation was predicted correctly - generate T/F matrix.

# Generate 1st - 2nd cosine score matrix

# Take scores for each protein and sum into families.

# Calculate T/F rates across the distribution for each family.

# Generate a confidence matrix. y-axis = families. x-axis = confidence scores at each point 0-1 in .01 increments.



