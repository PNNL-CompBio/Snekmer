"""features: Output feature processing for the Kmer pipeline.

author: @christinehc
"""

# imports
import os.path
from collections import Counter
from itertools import repeat
from multiprocessing import Pool
from kmerfeatures.transform import vectorize_string


# functions
def output_features(save_to, fformat, feature_sets=None,
                    labels=None, mode="w", **kwargs):
    """Generate output features based on input fasta file.

    Parameters
    ----------
    save_to : str
        File path to output file (default: None).
    fformat : str
        File format; select one option from the following list:
            ('gist', 'sieve', 'matrix', 'both')
    feature_sets : type
        Description of parameter `feature_sets`.
    labels : list
        Description of parameter `labels`.
    mode : str
        File write mode (default: 'w').
            For details, see documentation for Python built-in
            function open().

    Returns
    -------
    None
        Feature output is written to the `save_to` file path.
        No output is returned.

    """
    # # ensure valid format is specified
    # if fformat not in ['gist', 'sieve', 'matrix', 'both']:
    #     raise ValueError(
    #         ('File format must be one of the following:'
    #          ' "gist", "sieve", "both", or "matrix".')
    #         )

    # update any gist files
    if fformat in ["gist", "both"]:
        output_gist(save_to, labels=labels,
                    feature_sets=feature_sets, mode=mode, **kwargs)

    # update any sieve files
    if fformat in ["sieve", "both"]:
        output_sieve(save_to, feature_sets=feature_sets, mode=mode, **kwargs)

    # update matrix files
    if fformat == "matrix":
        output_matrix(save_to, labels=labels,
                      feature_sets=feature_sets, mode=mode)


def output_gist(filename, labels=None, feature_sets=None, mode='w'):
    """Output features in gist format.

    Parameters
    ----------
    filename : str
        /path/to/file.ext
    labels : list
        (default: None)
    feature_sets : type
        Description of parameter `feature_sets`.
    mode : str
        (default: 'w')

    Returns
    -------
    None
        Feature output is written to files; no output is returned.

    """
    train_out = f"{filename.split('.')[0]}.train"
    class_out = f"{filename.split('.')[0]}.class"

    # update train
    with open(train_out, mode) as tf, open(class_out, mode) as cf:
        tf.write("corner")
        cf.write("corner\tclass\n")

        if labels:
            for label in labels:
                tf.write("\t%s" % label)
        else:
            for i in range(len(feature_sets[0]) - 1):
                tf.write("\tlabel%d" % i)
        for features in feature_sets:
            output_gist_features(tf, features, mode)
            output_gist_class(tf, features, mode)


def output_sieve(filename, feature_sets=None, mode='w', **kwargs):
    """Output features in sieve format.

    Parameters
    ----------
    filename : str
        /path/to/file.ext
    labels : list
        (default: None)
    feature_sets : type
        Description of parameter `feature_sets`.
    mode : str
        (default: 'w')

    Returns
    -------
    None
        Feature output is written to file; no output is returned.

    """

    def output_sieve_features(features, filename, example_index=None):
        """Write features to SIEVE output file.

        Parameters
        ----------
        features : type
            List of specified features.
        filename : str
            File handle for naming of output files.
        example_index : dict
            Description of parameter `example_index` (default: None).

        Returns
        -------
        None
            Feature output is written to file; no output is returned.

        """
        # parse first item in feature list as feature ID
        fid = features[0]
        value = example_index.get(fid, 0.0)

        with open(filename, "a") as f:
            f.write("pattern\t%s\t%d\n" % (fid, len(features) - 1))
            f.write("\tinput\t%s" % fid)
            for ft in features[1:]:
                f.write("\t%s" % ft)
            f.write("\n")
            f.write("\toutput\t%s\t%d\n" % (fid, value))
            f.flush()

    # pattern_out = f"{filename}.pattern"
    with open(filename, mode) as f:
        for features in feature_sets:
            output_sieve_features(features, f, **kwargs)


def output_matrix(filename, labels=False, feature_sets=False, mode='w'):
    """Output features in matrix format.

    Parameters
    ----------
    filename : str
        /path/to/file.ext
    labels : list
        (default: False)
    feature_sets : type
        Description of parameter `feature_sets`.
    mode : str
        (default: 'w')

    Returns
    -------
    None
        Feature output is written to file; no output is returned.

    """
    with open(filename, mode) as f:
        if labels:
            f.write("%s" % labels[0])
            for label in labels[1:]:
                f.write("\t%s" % label)
            f.write("\n")
            f.flush()

    if feature_sets:
        for features in feature_sets:
            output_gist_features(filename, features, mode)


def output_gist_features(filename, features, mode='w'):
    """Write features to gist output file.

    Parameters
    ----------
    filename : type
        Description of parameter `file`.
    features : type
        Description of parameter `features`.

    Returns
    -------
    None
        Feature output is written to file; no output is returned.

    """
    with open(filename, mode) as f:
        f.write("%s" % features[0])
        for ft in features[1:]:
            f.write("\t%s" % ft)
        f.write("\n")
    # filename.flush()


def output_gist_class(filename, features, example_index=None, mode='w'):
    """Write gist class to specified output file.

    Parameters
    ----------
    filename : str
        /path/to/file.ext
    features : list
        List of specified features.
    example_index : dict
        Description of parameter `example_index` (default: None).

    Returns
    -------
    None
        Feature output is written to file; no output is returned.

    """
    pattern = f"{filename.split('.')[0]}.pattern"
    with open(pattern, mode) as f:
        fid = features[0]
        value = example_index.get(fid, -1)
        f.write("%s\t%d\n" % (features[0], value))
        # filename.flush()


def define_feature_space(sequence_dict, k, map_function=None, start=None,
                         end=None, min_rep_thresh=2, verbose=False,
                         log_file=False, processes=2):
    """Create a feature dictionary defined from parameters.

    Parameters
    ----------
    sequence_dict : dict
        Sequence dictionary {ID: sequence}.
    k : int
        Sequence length k of the kmer.
    map_function : str
        Name of the map function (e.g. "reduced_alphabet_0")
        (default: None; applies no mapping).
    start : int
        Start index of the sequence (for sequence slicing).
    end : int
        End index of the sequence (for sequence slicing).
    min_rep_thresh : type
        Description of parameter `min_rep_thresh`.
    verbose : bool
        If True, enables verbose output (default: False).
    log_file : str
        /path/to/log.file for verbose outputs (default: False)
        If False, pipes verbose outputs to console instead.
    processes : int
        (for multithreading) Number of partitions to create

    Returns
    -------
    dict
        Filtered feature dictionary {ID: sequence}.

    """
    # use multiprocessing to parallelize kmerization step
    feature_dict = {}
    with Pool(processes) as pool:
        feature_dict = pool.starmap(
            vectorize_string, zip(sequence_dict.values(),
                                  repeat(k),                # k
                                  repeat(start),            # start
                                  repeat(end),              # end
                                  repeat(map_function),     # map_function
                                  repeat(feature_dict),     # feature_dict
                                  repeat(None),             # filter_list
                                  repeat(None),             # exclude
                                  repeat(True),             # return_dict
                                  repeat(False),            # verbose
                                  repeat(False))            # log_file
            )
    # combine dictionaries and sum any values with common keys
    feature_dict = dict(sum(map(Counter, feature_dict), Counter()))

    # if this is between 0 and 1 then it's a percentage
    if 0 < min_rep_thresh < 1:
        min_rep_thresh = len(feature_dict.keys()) * min_rep_thresh

    # filter out all those below the min_rep_thresh
    if min_rep_thresh:
        filter_dict = {}
        for key in feature_dict.keys():
            if feature_dict[key] >= min_rep_thresh:
                filter_dict[key] = feature_dict[key]
    else:
        filter_dict = feature_dict

    if verbose and log_file:
        with open(log_file, 'w') as f:
            f.write(
                ("Feature space: {0} kmers with more than "
                 "{1} representation in {2} sequences").format(
                     len(filter_dict),
                     min_rep_thresh,
                     len(sequence_dict)
                     )
                 )
    elif verbose and not log_file:
        print(
            ("Feature space: {0} kmers with more than "
             "{1} representation in {2} sequences").format(
                 len(filter_dict),
                 min_rep_thresh,
                 len(sequence_dict)
                 )
             )

    filter_list = filter_dict.keys()
    if len(filter_list) == 0:
        raise ValueError(
            ("Prefiltered feature space cannot be empty.")
            )

    return filter_dict
