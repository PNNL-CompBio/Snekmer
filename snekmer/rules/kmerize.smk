"""kmerize.smk: Module for kmer vector generation.

author: @christinehc
"""

# include unzipping module
include: "process_input.smk"


# built-in imports
import gzip
import json
from datetime import datetime
from glob import glob
from itertools import product, repeat
from os.path import basename, join

# external libraries
import snekmer as skm

# get input files
input_files = glob(join("input", "*"))
unzipped = [
    fa.rstrip(".gz")
    for fa, ext in product(input_files, ["fasta"])
    if fa.rstrip(".gz").endswith(f".{ext}")
]
zipped = [fa for fa in input_files if fa.endswith(".gz")]
UZS = [skm.utils.split_file_ext(f)[0] for f in zipped]
FAS = [skm.utils.split_file_ext(f)[0] for f in unzipped]


# rule all:
#     input:
#         expand(join("input", '{uz}'), uz=UZS),  # require unzipping
#         expand(join("output", "features", "{nb}", "{fa}.json"), nb=NON_BGS, fa=FAS)


rule generate:
    input:
        params=join("output", "processed", "{nb}.json"),
    output:
        labels=join("output", "labels", "{nb}.txt"),
    log:
        join("output", "labels", "log", "{nb}.log"),
    run:
        start_time = datetime.now()

        # read processed features
        with open(input.params, "r") as f:
            params = json.load(f)

        # generate labels only
        labels = skm.transform.generate_labels(
            config["k"],
            alphabet=params["alphabet"],
            filter_list=params["filter_list"],
        )
        if config["output"]["format"] == "simple":
            skm.features.output_features(output.labels, "matrix", labels=labels)

        # record script runtime
        skm.utils.log_runtime(log[0], start_time)



rule vectorize:
    input:
        kmers=join("output", "labels", "{nb}.txt"),
        params=join("output", "processed", "{nb}.json"),
        fastas=join("input", "{nb}.fasta"),
    log:
        join("output", "features", "log", "{nb}.log"),
    output:
        files=expand(join("output", "features", "{{nb}}", "{fa}.json.gz"), fa=FAS),
    run:
        start_time = datetime.now()

        # get kmers for this particular set of sequences
        kmers = skm.io.read_output_kmers(input.kmers)

        # read processed features
        with open(input.params, "r") as f:
            params = json.load(f)

        # sort i/o lists to match wildcard order
        fastas = sorted(input.fastas)
        outfiles = sorted(output.files)

        # revectorize based on full kmer list
        for i, fa in enumerate(fastas):
            results = {"seq_id": [], "vector": []}
            seq_list, id_list = skm.io.read_fasta(fa)
            for seq, sid in zip(seq_list, id_list):
                results["seq_id"] += [sid]
                results["vector"] += [
                    skm.transform.vectorize_string(
                        seq,
                        config["k"],
                        params["alphabet"],
                        start=config["start"],
                        end=config["end"],
                        filter_list=kmers,  # params['filter_list'],
                        verbose=False,  # way too noisy for batch process
                        log_file=log[0],
                    )
                ]

            # with open(outfiles[i], 'w') as f:
            #     json.dump(results, f)

            with gzip.open(outfiles[i], "wt", encoding="ascii") as zipfile:
                json.dump(results, zipfile)

        # record script runtime
        skm.utils.log_runtime(log[0], start_time)

rule vectorize_full:
    input:
        kmers=join("output", "labels", "full", "{nb}.txt"),
        params=join("output", "processed", "{nb}.json"),
        fasta=join("input", "{nb}.fasta"),
    log:
        join("output", "features", "log", "{nb}.log"),
    output:
        file=join("output", "features", "full", "{nb}.json.gz"),
    run:
        start_time = datetime.now()

        # get kmers for this particular set of sequences
        kmers = skm.io.read_output_kmers(input.kmers)

        # read processed features
        with open(input.params, "r") as f:
            params = json.load(f)

        # revectorize based on full kmer list
        # for i, fa in enumerate(fastas):
        results = {"seq_id": [], "vector": []}
        seq_list, id_list = skm.io.read_fasta(input.fasta)
        for seq, sid in zip(seq_list, id_list):
            results["seq_id"] += [sid]
            results["vector"] += [
                skm.transform.vectorize_string(
                    seq,
                    config["k"],
                    params["alphabet"],
                    start=config["start"],
                    end=config["end"],
                    filter_list=kmers,  # params['filter_list'],
                    verbose=False,  # way too noisy for batch process
                    log_file=log[0],
                )
            ]

        with gzip.open(output.file, "wt", encoding="ascii") as zipfile:
            json.dump(results, zipfile)

        # record script runtime
        skm.utils.log_runtime(log[0], start_time)

rule vectorize_search:
    input:
        #kmers=join("output", "labels", "{nb}.txt"),
        #params=join("output", "processed", "{nb}.json"),
        fastas=join("input", "{nb}.fasta"),
    log:
        join("output", "features", "log", "{nb}.log"),
    output:
        files=expand(join("output", "features", "{{nb}}", "{fa}.json.gz"), fa=FAS),
    run:
        start_time = datetime.now()

        # get kmers for this search
        kmers = []
        with open(config["input"]["feature_set"], "r") as f:
            kmers = skm.io.read_output_kmers(config["input"]["feature_set"])

        # old
        #kmers = skm.io.read_output_kmers(input.kmers)

        # read processed features
        #with open(input.params, "r") as f:
        #    params = json.load(f)

        # sort i/o lists to match wildcard order
        fastas = sorted(input.fastas)
        # print(fastas)
        outfiles = sorted(output.files)

        # revectorize based on full kmer list
        for i, fa in enumerate(fastas):
            results = {"seq_id": [], "vector": []}
            seq_list, id_list = skm.io.read_fasta(fa)
            for seq, sid in zip(seq_list, id_list):
                results["seq_id"] += [sid]
                results["vector"] += [
                    skm.transform.vectorize_string(
                        seq,
                        config["k"],
                        config["alphabet"],
                        start=config["start"],
                        end=config["end"],
                        filter_list=kmers,  # params['filter_list'],
                        verbose=False,  # way too noisy for batch process
                        log_file=log[0],
                    )
                ]

            # with open(outfiles[i], 'w') as f:
            #     json.dump(results, f)

            with gzip.open(outfiles[i], "wt", encoding="ascii") as zipfile:
                json.dump(results, zipfile)

        # record script runtime
        skm.utils.log_runtime(log[0], start_time)
