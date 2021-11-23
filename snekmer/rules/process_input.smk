"""process_input.smk: Module for input file handling.

author: @christinehc
"""

# imports
from datetime import datetime
from os import makedirs
from os.path import dirname, exists, join
from pandas import DataFrame
import snekmer as skm

# define rules
rule unzip:
    input:
        join("input", "{uz}.gz")
        # lambda wildcards: join("input", f"{wildcards.uz}.{uz_map[wildcards.uz]}.gz")
        # lambda wildcards: unzip(join("input", ))
    output:
        join("input", "{uz}")
    params:
        outdir=join("input", 'zipped')
    shell:
        "mkdir {params.outdir} && cp {input} {params.outdir} && gunzip -c {input} > {output}"
    # run:
        # if wildcards.sample.endswith('.fastq'):
        #     shell("echo gzip {input}")
        #     shell("echo mv {input}.gz {params.outdir}")
        # else:
        #     shell("mv {input} {params.outdir}")


# rule perform_kmer_walk:
#     input:
#         fasta=get_fasta_files
#     output:
#         # need to fix code to properly generate an output...
#     run:
#         skm.walk.kmer_walk(input.fasta)


# read and process parameters from config
rule preprocess:
    input:
        fasta=lambda wildcards: join("input", f"{wildcards.nb}.{fa_map[wildcards.nb]}"),
    output:
        data=join("output", "processed", "{nb}.json"),
        desc=join("output", "processed", "{nb}_description.csv")
    log:
        join("output", "processed", "log", "{nb}.log")
    run:
        # log step initialization
        start_time = datetime.now()

        # read fasta file
        seq_list, id_list = skm.io.read_fasta(input.fasta)

        # if random alphabet specified, implement randomization
        if config['randomize_alphabet']:
            rand_alphabet = skm.transform.randomize_alphabet(config['input']['alphabet'])
            alphabet = [residues, map_name, rand_alphabet]
        else:
            alphabet = config['alphabet']
            if alphabet == "None":
                alphabet = None

        # if no feature set is specified, define feature space
        if not config['input']['feature_set']:
            # prefilter fasta to cut down on the size of feature set
            filter_dict = skm.features.define_feature_space(
                {k: v for k, v in zip(id_list, seq_list)},
                config['k'],
                alphabet=alphabet,
                start=config['start'],
                end=config['end'],
                min_rep_thresh=min_rep_thresh,
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
                filter_list = skm.io.read_output_kmers(config['input']['feature_set'])

        # optional indexfile with IDs of good feature output examples
        if config['input']['example_index_file']:
            example_index = skm.io.read_example_index(
                config['input']['example_index_file']
            )
        else:
            example_index = {}

        # loop thru seqs, apply input params to preprocess seq list
        seen = []  # filter duplicates
        save_data = dict()

        # define recursive and nonrecursive saving patterns for params
        recursive = ['sequences', 'ids', 'residues']
        nonrecursive = ['alphabet', 'k', 'example_index', 'filter_list']
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
                scid_list, scramble_list, example_index = skm.transform.scramble_sequence(
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
                sids_n, seqs_n = skm.transform.make_n_terminal_fusions(
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
            [alphabet, config['k'], example_index, filter_list],
            nonrecursive
        ):
            save_data[label] = dset

        # save all parameters into json
        with open(output.data, 'w') as f:
            json.dump(save_data, f)

        # read and save fasta descriptions into dataframe
        try:
            desc = skm.utils.parse_fasta_description(input.fasta)
            desc.to_csv(output.desc)
        except AttributeError:  # if no description exists > empty df
            DataFrame([]).to_csv(output.desc)

        # record script runtime
        skm.utils.log_runtime(log[0], start_time)
