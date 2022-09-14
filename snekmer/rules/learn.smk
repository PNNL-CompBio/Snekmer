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

# collect all fasta-like files, unzipped filenames, and basenames
input_dir = "input" if (("input_dir" not in config) or (str(config["input_dir"]) == "None")) else config["input_dir"]
input_files = glob(join(input_dir, "*"))
zipped = [fa for fa in input_files if fa.endswith(".gz")]
unzipped = [
    fa.rstrip(".gz")
    for fa, ext in product(input_files, config["input_file_exts"])
    if fa.rstrip(".gz").endswith(f".{ext}")
]

annot_files = glob(join("annotations", "*.ann"))

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

# define output files to be created by snekmer
rule all:
    input:
        #expand(join(input_dir, "{uz}"), uz=UZS),  # require unzipping
        join("output", "learn", "kmer-associations.csv"),  # require learning


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
        data=expand(join("output", "vector", "{fa}.npz"), fa=NON_BGS),
        annotation=expand("{an}", an=annot_files),
    output:
        table=join("output", "learn", "kmer-associations.csv"),
        summary=join("output", "learn", "kmer-summary.csv"),
    log:
        join(out_dir, "learn", "log", "learn.log"),
    run:
        # log script start time
        start_time = datetime.now()
        with open(log[0], "a") as f:
            f.write(f"start time:\t{start_time}\n")

        #load in the kmerbasis and vectors
        for f in input.data:
            df, kmerlist = skm.io.load_npz(f)
            seqids = df["sequence_id"]
            vecs = skm.utils.to_feature_matrix(df["sequence_vector"].values)
            
        def count_kmer(sequence, kmer):
            count = start = 0
            while True:
                start = sequence.find(kmer, start) + 1
                if start > 0:
                    count+=1
                else:
                    return count

        # vecs is now an numpy array with columns corresponding to kmerlist
        #      each row is a protein with names given by seqids
        # seqids: >tr|A0A0C2DJP9|A0A0C2DJP9_9STAP Lipid kinase OS=Salinicoccus roseus OX=45670 GN=SN16_09640 PE=3 SV=1

        # read in the annotation files which should be in the format
        #    column 1: sequence id
        #    column 2: annotation
        # annotations: A0A0C2DJP9	IPR016064

        # The plan is to build a matrix of kmer-annotation associations
        #    by calculating the probability that a kmer is associated with
        #    each annotation (label). Options should give the ability to
        #    filter out low information associations from the output
        # And then output an output/learn/kmer-associations.csv file (or some
        #    format file)


        annots = list()
        for f in input.annotation:
            annots.append(pd.read_table(f))
        annotations = pd.concat(annots)

        klist_length = len(kmerlist)
        seq_kmer_dict = {}
        counter = 0
        for seq in seqids:
            seq_kmer_dict[seq] = [None]*klist_length
            seq_kmer_dict[seq].insert(0,(df["sequence"][counter]))
            counter +=1

        kmer_totals = {}
        for item in kmerlist:
            kmer_totals[item] = 0

        for seq in seq_kmer_dict:
            count=1
            for kmer in kmerlist:
                kmers_in_seq = count_kmer(str(seq_kmer_dict[seq][0]),str(kmer))
                seq_kmer_dict[seq][count] = kmers_in_seq
                kmer_totals[kmer] += kmers_in_seq
                count+=1
                
        for seq in seq_kmer_dict:
            count=1
            for kmer in kmerlist:
                seq_kmer_dict[seq][count] = (seq_kmer_dict[seq][count]/kmer_totals[kmer])
                count+=1


        seq_kmer_dict = {seq: seq_kmer_dict[seq][1:] for seq in seq_kmer_dict}

        #probability vector will be saved as kmer-associations.csv
        probability_vector = []
        for value in seq_kmer_dict.values():
            probability_vector.append(value)



        # Filter Methods - optional
        # Filter by value = 0
        # Filter by percentile = 1
        # Filter by rank = 2

        filter_method = 2

        filtered_pv = copy.deepcopy(probability_vector)

        #keep kmers greater than N probability value
        if filter_method == 0:
            filter_by_value = .01
            for i,row in enumerate(probability_vector):
                for j,column in enumerate(row):
                    if column < filter_by_value:
                        filtered_pv[i][j] = 0

        #Filter out lower N percent of kmers
        if filter_method == 1:
            filter_by_percent = 90
            pv_np = np.array(probability_vector).astype("float")
            pv_np[pv_np==0] =np.nan
            for i,row in enumerate(probability_vector):
                filt = np.nanpercentile(pv_np[i],filter_by_percent)
                for j,column in enumerate(row):
                    if column < filt:
                        filtered_pv[i][j] = 0
        #Filter kmers not in top N 
        if filter_method == 2:
            filter_by_rank = 10
            pv_np = np.array(probability_vector)
            ranked_pv = copy.deepcopy(probability_vector)
            for i,row in enumerate(probability_vector):
                ranked_pv[i] = rankdata(row, method="ordinal")
            for i,row in enumerate(probability_vector):
                for j,column in enumerate(row):
                    if ranked_pv[i][j] <= (len(row)-filter_by_rank):
                        filtered_pv[i][j] = 0


        #summary vector for kmer-summary.csv
        summary_vect = []
        for i,row in enumerate(filtered_pv):
            row_summary = []
            row_summary.append(list(seq_kmer_dict.keys())[i])
            for j,value in enumerate(row):
                if value != 0:
                    row_summary.append(kmerlist[j])
                    row_summary.append(value)
            summary_vect.append(row_summary)
                

        with open("output/learn/kmer-associations.csv", "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerows(probability_vector)
              
        # Summary Format based on filter: (SeqID: kmer,weight,kmer,weight,etc)
        with open("output/learn/kmer-summary.csv", "w", newline="") as fl:
            writer = csv.writer(fl)
            writer.writerows(summary_vect)

        skm.utils.log_runtime(log[0], start_time)

      #Validation - Due to rounding, these are not exactly equal to 1, but very close
        # sums_vector = [0]*klist_length
        # for i,value in enumerate(probability_vector):
        #     for index in range(0,klist_length):
        #         sums_vector[index] += value[index]
        # print(sums_vector)