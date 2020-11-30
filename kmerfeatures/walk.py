"""walk: Kmer walk module for Kmer pipeline.

author: @christinehc
"""
# imports
from kmerfeatures.utils import read_fasta
from kmerfeatures.transform import string_vectorize

# functions
def kmer_walk(fastafile, max_k=20, seq_dump=False):
    """Short summary.

    Parameters
    ----------
    fastafile : type
        Description of parameter `fastafile`.
    max_k : type
        Description of parameter `max_k`.
    seq_dump : type
        Description of parameter `seq_dump`.

    Returns
    -------
    type
        seems to return print statements
        can modify this to output to some log file.

    """
    # next read in sequences from the fasta file
    seq_list, id_list = read_fasta(fastafile, include_map=False)

    # eliminate the bottom n% of kmers by representation
    min_thresh = 10

    exclusion_list = None
    for kmer in range(1, max_k):
        # build a feature dict for this kmer and these sequences
        feature_dict = {}

        # short circuit the exclusion list (very slow) but allow it to count
        exclusion_list = []
        for i in range(len(seq_list)):
            seq = seq_list[i]
            seq_id = id_list[i]

            feature_dict = string_vectorize(sequence=seq, kmer=kmer,
                                            map_function="reduced_alphabet_0",
                                            feature_dict=feature_dict,
                                            exclusion_list=exclusion_list,
                                            return_dict=True)

        exclusion_list = []
        for key in feature_dict.keys():
            if feature_dict[key] < min_thresh:
                exclusion_list.append(key)
        print("Kmer %d, number of remaining features %d total, number of remaining features occuring more than %d times %d, of %g possible, %g%%" %
              (kmer, len(feature_dict.keys()), min_thresh, len(feature_dict.keys())-len(exclusion_list), 20**kmer, (len(feature_dict.keys())-len(exclusion_list))/20**kmer))

    if seq_dump:
        for seq in feature_dict.keys():
            if feature_dict[seq] > min_thresh:
                print(seq)
