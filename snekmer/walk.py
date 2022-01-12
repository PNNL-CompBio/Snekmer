"""walk: Snekmer walk module.

author: @christinehc
"""
# imports
from .utils import read_fasta
from .transform import vectorize_string


# functions
def kmer_walk(fasta, max_k=20, seq_dump=False, min_thresh=10,
              map_function='reduced_alphabet_0'):
    """Short summary.

    Parameters
    ----------
    fasta : str
        path/to/fasta_file.fasta
    max_k : int
        Maximum kmer value (default: 20).
    seq_dump : bool
        Description of parameter `seq_dump`.
    min_thresh : float
        Representation % threshold, below which kmers are eliminated
        (default: 10).
    map_function : str
        Reduced alphabet name in the format 'reduced_alphabet_n',
        where 0 <= n <= 4 (default: 'reduced_alphabet_0').
        (See __init__.py or SIEVEInit for map descriptions.)

    Returns
    -------
    type
        seems to return print statements
        can modify this to output to some log file.

    """
    # read sequences from fasta file
    seq_list, _ = read_fasta(fasta)

    exclusion_list = None
    for kmer in range(1, max_k):
        # build feature dict for this kmer and these sequences
        feature_dict = {}
        exclusion_list = []

        # short-circuit the exclusion list (very slow) but allow it to count
        for seq in seq_list:
            feature_dict = vectorize_string(seq, kmer,
                                            map_function=map_function,
                                            feature_dict=feature_dict,
                                            exclude=exclusion_list,
                                            return_dict=True)

        exclusion_list = []
        for key in feature_dict.keys():
            if feature_dict[key] < min_thresh:
                exclusion_list.append(key)

        print(f"Kmer {kmer}\nnumber of remaining features:"
              f" {len(feature_dict.keys())} total\n"
              "number of remaining features occuring"
              f" more than {min_thresh} times"
              f" {len(feature_dict) - len(exclusion_list)},"
              f" of {20**kmer} possible,"
              f" {(len(feature_dict) - len(exclusion_list)) / 20**kmer}"
              )


    if seq_dump:
        for key in feature_dict.keys():
            if feature_dict[key] >= min_thresh:
                print(seq)
