# built-in imports
import glob
import json
from datetime import datetime
from itertools import (product, repeat)
from multiprocessing import Pool
from os import makedirs
from os.path import (basename, dirname, exists, join, splitext)

# external libraries
from kmerfeatures import (alphabet, features, score, transform, utils, walk)
import numpy as np
from pandas import DataFrame
from Bio import SeqIO


# collect all fasta-like files, unzipped filenames, and basenames
exts = ['fasta', 'fna', 'faa', 'fa']
input_files = glob.glob(join(config['input']['fasta_dir'], "*"))
compressed = [fa for fa in input_files if fa.endswith('.gz')]
uncompressed = [fa.rstrip('.gz') for fa, ext in product(input_files, exts) if fa.rstrip('.gz').endswith(f".{ext}")]

unzipped_ext_map = {splitext(splitext(basename(f))[0])[0]: splitext(splitext(basename(f))[0])[1].lstrip('.') for f in compressed}
fasta_ext_map = {splitext(basename(f))[0]: splitext(basename(f))[1].lstrip('.') for f in uncompressed}
UZS = list(unzipped_ext_map.keys())
FAS = list(fasta_ext_map.keys())
print(fasta_ext_map)
print(FAS)

# define output directory (helpful for multiple runs)
out_dir = join(config['output']['save_dir'],
               config['map_function'],
               f"k-{config['k']:02}")

# validity check
if config['map_function'] not in alphabet.ALPHABETS.keys():
    raise ValueError("Invalid alphabet specified; alphabet must be a"
                     " string in the form 'reduced_alphabet_n' where"
                     " n is an integer between 0 and"
                     f" {len(alphabet.ALPHABETS)}.")


# define workflow
rule all:
    input:
        # expand(join(config['input']['fasta_dir'], '{uz}.{ext}'), uz=UNZIPPED.keys(), ext=UNZIPPED.values()),
        expand(join(out_dir, "features", "{fa}", "{fa2}.json"), fa=FAS, fa2=FAS),
        # expand(join(out_dir, "features", "{fa}.json"), fa=FAS.keys()),
        expand(join(out_dir, "score", "{fa}.json"), fa=FAS)


# [in-progress] kmer walk
if config['walk']:
    rule perform_kmer_walk:
        input:
            fasta=get_fasta_files
        output:
            # need to fix code to properly generate an output...
        run:
            walk.kmer_walk(input.fasta)


# # if any files are gzip compressed, unzip them
# if len(UZS) > 0:
#     rule unzip:
#         input:
#             # lambda wildcards: '{0}_hmm'.format(dict[wildcards.run])
#             lambda wildcards: join(config['input']['fasta_dir'], f"{wildcards.uz}.{unzipped_ext_map[wildcards.uz]}.gz")
#         output:
#             join(config['input']['fasta_dir'], "{uz}.{uzext}")
#         params:
#             outdir=join(config['input']['fasta_dir'], 'compressed')
#         shell:
#             "mkdir {params.outdir} && cp {input} {params.outdir} && gunzip -c {input} > {output}"


# read and process parameters from config
rule preprocess:
    input:
        fasta=lambda wildcards: join(config['input']['fasta_dir'], f"{wildcards.fa}.{fasta_ext_map[wildcards.fa]}"),
    output:
        data=join(out_dir, "processed", "{fa}.json"),
        desc=join(out_dir, "processed", "{fa}_description.json")
    log:
        join(out_dir, "processed", "log", "{fa}.log")
    run:
        # log step initialization
        start_time = datetime.now()

        # read fasta file
        seq_list, id_list = utils.read_fasta(input.fasta)

        # if random alphabet specified, implement randomization
        if config['randomize_alphabet']:
            rand_alphabet = transform.randomize_alphabet(config['input']['map_function'])
            map_function = [residues, map_name, rand_alphabet]
        else:
            map_function = config['map_function']

        # if no feature set is specified, define feature space
        if not config['input']['feature_set']:
            # prefilter fasta to cut down on the size of feature set
            filter_dict = features.define_feature_space(
                {k: v for k, v in zip(id_list, seq_list)},
                config['k'],
                map_function=map_function,
                start=config['start'],
                end=config['end'],
                min_rep_thresh=config['min_rep_thresh'],
                verbose=config['output']['verbose'],
                log_file=log[0],
                processes=config['processes']
                )
            filter_list = list(filter_dict.keys())
            assert len(filter_list) > 0, "Invalid feature space; terminating."
        else:
            # read in list of ids to use from file; NO FORMAT CHECK
            filter_list = []
            with open(config['input']['feature_set'], "r") as f:
                for line in f.readlines():
                    filter_list.append(line.split()[0])

        # optional indexfile with IDs of good feature output examples
        if config['input']['example_index_file']:
            example_index = utils.read_example_index(
                config['input']['example_index_file']
                )
        else:
            example_index = {}

        # loop thru seqs, apply input params to preprocess seq list
        seen = []  # filter duplicates
        save_data = dict()

        # define recursive and nonrecursive saving patterns for params
        recursive = ['sequences', 'ids', 'residues']
        nonrecursive = ['map_function', 'k', 'example_index', 'filter_list']
        all_dsets = recursive + nonrecursive

        for i in range(len(seq_list)):
            seq = seq_list[i]
            sid = id_list[i]

            # ignore duplicate ids
            if config['output']['filter_duplicates'] and sid in seen:
                continue
            seen.append(sid)

            seqs = [seq]
            sids = [sid]

            # shuffle the N-terminal sequence n times
            if config['output']['shuffle_n']:
                example_index[id] = 1.0
                scid_list, scramble_list, example_index = transform.scramble_sequence(
                    sid, seq[:30], n=config['output']['shuffle_n'],
                    example_index=example_index
                    )
                seqs += scramble_list
                sids += scid_list

                # include shuffled sequences in output
                if config['output']['shuffle_sequences']:
                    filename = join(out_dir, 'shuffled',
                                    wildcards.fa, "%s_shuffled.fasta" % sid)
                    if not exists(dirname(filename)):
                        makedirs(dirname(filename))
                    with open(filename, "w") as f:
                        for i in range(len(sids)):
                            f.write(">%s\n%s\n" % (sids[i], seqs[i]))

            # run SIEVE on the wt and each shuffled sequence
            if config['output']['n_terminal_file']:
                sids_n, seqs_n = transform.make_n_terminal_fusions(
                    sid, config['output']['n_terminal_file']
                    )
                seqs += seqs_n
                sids += sids_n
            residues = None
            if config['nucleotide']:
                residues = "ACGT"

            # populate dictionary for json save file
            to_save = [seqs, sids, residues]
            save_label = recursive
            for dset, label in zip(to_save, save_label):
                if label in save_data.keys() and save_data[label] is not None:
                    save_data[label] = save_data[label] + dset
                else:
                    save_data[label] = dset

        # save variables not generated in loop
        for dset, label in zip(
            [map_function, config['k'], example_index, filter_list],
            nonrecursive
        ):
            save_data[label] = dset

        # save all parameters into json
        with open(output.data, 'w') as f:
            json.dump(save_data, f)

        # read and save fasta descriptions into dataframe
        try:
            desc = utils.parse_fasta_description(input.fasta)
            desc.to_json(output.desc)
        except AttributeError:  # if no description exists > empty df
            DataFrame([]).to_json(output.desc)

        # record script runtime
        end_time = datetime.now()
        with open(log[0], 'a') as f:
            f.write(f"start time:\t{start_time}\n")
            f.write(f"end time:\t{end_time}\n")
            f.write(f"total time:\t{utils.format_timedelta(end_time - start_time)}")


rule generate_kmer_labels:
    input:
        params=rules.preprocess.output.data
    output:
        labels=join(out_dir, "labels", "{fa}.txt")
    log:
        join(out_dir, "labels", "log", "{fa}.log")
    run:
        start_time = datetime.now()

        # read processed features
        with open(input.params, 'r') as f:
            params = json.load(f)

        # generate labels only
        labels = transform.generate_labels(
            config['k'],
            map_function=params['map_function'],
            residues=params['residues'],
            filter_list=params['filter_list']
            )
        if config['output']['format'] == "simple":
            features.output_features(
                output.labels, "matrix", labels=labels
                )

        # record script runtime
        end_time = datetime.now()
        with open(log[0], 'a') as f:
            f.write(f"start time:\t{start_time}\n")
            f.write(f"end time:\t{end_time}\n")
            f.write(f"total time:\t{utils.format_timedelta(end_time - start_time)}")


# def get_fasta_chunks(wildcards):
#     files = glob.glob(
#         join(out_dir, "processed", f"{wildcards.fa}", f"{wildcards.fa}_*.fasta")
#         )
#     return files


# rule generate_kmer_features:
#     input:
#         params=rules.preprocess.output.data,
#         labels=rules.generate_kmer_labels.output.labels
#     output:
#         features=join(out_dir, "features", "{fa}.txt")
#     log:
#         join(out_dir, "features", "log", "{fa}.log")
#     run:
#         start_time = datetime.now()
#
#         # read processed features
#         with open(input.params, 'r') as f:
#             params = json.load(f)
#
#         # apply user-specified save name, if it exists
#         # if config['output']['filename'] is None:
#         #     output_file = wildcards.fa
#
#         # generate features for each sequence and output features
#         first = True
#         for i in range(len(params['sequences'])):
#             seq = params['sequences'][i]
#             seq_id = params['ids'][i]
#
#             # labels = []
#
#             if config['output']['verbose']:
#                 with open(log[0], 'a') as f:
#                     f.write(f"Constructing features for sequence {seq_id}\n")
#
#             feature_list = [seq_id]
#
#             feature_list += transform.vectorize_string(
#                 seq,
#                 k=config['k'],
#                 start=config['start'],
#                 end=config['end'],
#                 map_function=params['map_function'],
#                 filter_list=params['filter_list'],  # utils.get_output_kmers(input.labels),
#                 verbose=config['output']['verbose'],
#                 log_file=log[0]
#                 )
#
#             # record labels for first sequence only
#             if first:
#                 labels = utils.get_output_kmers(input.labels)
#                 # labels += transform.generate_labels(
#                 #     config['k'],
#                 #     map_function=params['map_function'],
#                 #     residues=params['residues'],
#                 #     filter_list=params['filter_list']
#                 #     )
#                 if config['output']['format'] == "simple":
#                     features.output_features(output.features,
#                                                  "matrix",
#                                                  labels=labels)
#
#             first = False
#
#             # output as we go (esp. good for very large input files)
#             if config['output']['format'] == "simple":
#                 features.output_features(
#                     output.features, "matrix", feature_sets=[features],
#                     mode="a"
#                     )
#
#             # output sieve patterns as we go to provide a record
#             if config['output']['format'] in ("sieve", "both"):
#                 features.output_features(
#                     output.features, "sieve", feature_sets=[features],
#                     mode="a", example_index=params['example_index']
#                     )
#
#             # only append features if not dumping into file
#             if config['output']['format'] != "simple":
#                 feature_sets.append(features)
#                 features.output_features(
#                     output.features,
#                     config['output']['format'],
#                     feature_sets=feature_sets,
#                     example_index=params['example_index'],
#                     labels=labels)
#
#         # record script runtime
#         end_time = datetime.now()
#         with open(log[0], 'a') as f:
#             f.write(f"start time:\t{start_time}\n")
#             f.write(f"end time:\t{end_time}\n")
#             f.write(f"total time:\t{utils.format_timedelta(end_time - start_time)}")


def get_expanded_fastas(wildcards):
    return expand(join(config['input']['fasta_dir'], f"{{fa2}}.{fasta_ext_map[fa2]}"),
                  fa2=FAS)


rule standardize_kmers:
    input:
        kmers=rules.generate_kmer_labels.output.labels,
        params=rules.preprocess.output.data,
        fastas=get_expanded_fastas
    log:
        join(out_dir, "features", "log", "{fa}.log")
    output:
        files=expand(join(out_dir, "features", "{{fa}}", "{fa2}.json"),
                     fa2=FAS)
    run:
        start_time = datetime.now()

        # get kmers for this particular set of sequences
        kmers = utils.get_output_kmers(input.kmers)

        # read processed features
        with open(input.params, 'r') as f:
            params = json.load(f)

        # sort i/o lists to match wildcard order
        fastas = sorted(input.fastas)
        outfiles = sorted(output.files)

        # revectorize based on full kmer list
        results = {'seq_id': [], 'vector': []}
        for i, fa in enumerate(fastas):
            seq_list, id_list = utils.read_fasta(fa)
            for seq, sid in zip(seq_list, id_list):
                results['seq_id'] += [sid]
                results['vector'] += [transform.vectorize_string(
                    seq,
                    k=config['k'],
                    start=config['start'],
                    end=config['end'],
                    map_function=params['map_function'],
                    filter_list=kmers,  # params['filter_list'],
                    verbose=False,  # way too noisy for batch process
                    log_file=log[0]
                )]

            with open(outfiles[i], 'w') as f:
                json.dump(results, f)

        # record script runtime
        end_time = datetime.now()
        with open(log[0], 'a') as f:
            f.write(f"start time:\t{start_time}\n")
            f.write(f"end time:\t{end_time}\n")
            f.write(f"total time:\t{utils.format_timedelta(end_time - start_time)}")


rule score_features:
    input:
        files=expand(join(out_dir,
                          "features", "{{fa}}", "{fa2}.json"),
                     fa2=FAS),
        fasta=join(config['input']['fasta_dir'], "{fa}.fasta")
    output:
        df=join(out_dir, "features", "{fa}.json"),
        scores=join(out_dir, "score", "{fa}.json")
    log:
        join(out_dir, "score", "log", "{fa}.log")
    run:
        start_time = datetime.now()
        with open(log[0], 'a') as f:
            f.write(f"start time:\t{start_time}\n")

        # parse all data
        label = config['score']['lname']
        data = utils.vecfiles_to_df(
            input.files, labels=config['score']['labels'], label_name=label
            )

        timepoint = datetime.now()
        with open(log[0], 'a') as f:
            f.write(f"vecfiles_to_df time:\t{utils.format_timedelta(timepoint - start_time)}\n")

        # parse family names and only add if some are valid
        families = [utils.get_family(fn) for fn in data['filename']]
        if any(families):
            label = 'family'
            data[label] = families

        # define feature matrix of kmer vectors
        feature_matrix = score.to_feature_matrix(data['vector'].values)

        # compute class probabilities
        labels = data[label].values
        class_probabilities = score.feature_class_probabilities(
            feature_matrix.T, labels, processes=config['processes']
            ).drop(columns=['presence'])

        new_timepoint = datetime.now()
        with open(log[0], 'a') as f:
            f.write(f"class_probabilities time:\t{utils.format_timedelta(new_timepoint - timepoint)}\n")

        # [IN PROGRESS] compute clusters
        clusters = score.cluster_feature_matrix(feature_matrix)
        data['cluster'] = clusters

        # save all files to respective outputs
        data.to_json(output.df)
        # np.save(output.npy, feature_matrix)
        class_probabilities.to_json(output.scores)

        # record script runtime
        end_time = datetime.now()
        with open(log[0], 'a') as f:
            f.write(f"total time:\t{utils.format_timedelta(end_time - start_time)}")
