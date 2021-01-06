"""walk: Kmer walk module for Kmer pipeline.

author: @christinehc
"""
# imports
from kmerfeatures.utils import read_fasta
from kmerfeatures.transform import vectorize_string

# functions
def kmer_walk(fastafile, max_k=20, seq_dump=False, min_thresh=10):
    """Short summary.

    Parameters
    ----------
    fastafile : type
        Description of parameter `fastafile`.
    max_k : type
        Description of parameter `max_k`.
    seq_dump : type
        Description of parameter `seq_dump`.
    min_thresh : float
        Representation % threshold, below which kmers are eliminated
        (default: 10).

    Returns
    -------
    type
        seems to return print statements
        can modify this to output to some log file.

    """
    # read sequences from fasta file
    seq_list, id_list = read_fasta(fastafile)

    exclusion_list = None
    for kmer in range(1, max_k):
        # build feature dict for this kmer and these sequences
        feature_dict = {}

        # short-circuit the exclusion list (very slow) but allow it to count
        exclusion_list = []
        for i in range(len(seq_list)):
            seq = seq_list[i]
            seq_id = id_list[i]

            feature_dict = vectorize_string(sequence=seq, kmer=kmer,
                                            map_function="reduced_alphabet_0",
                                            feature_dict=feature_dict,
                                            exclusion_list=exclusion_list,
                                            return_dict=True)

        exclusion_list = []
        for key in feature_dict.keys():
            if feature_dict[key] < min_thresh:
                exclusion_list.append(key)
        print("Kmer %d, number of remaining features %d total, number of remaining features occuring more than %d times %d, of %g possible, %g%%" %
              (kmer, len(feature_dict.keys()), min_thresh,
               len(feature_dict) - len(exclusion_list),
               20**kmer,
               (len(feature_dict) - len(exclusion_list)) / 20**kmer))

    if seq_dump:
        for key in feature_dict.keys():
            if feature_dict[key] >= min_thresh:
                print(seq)
