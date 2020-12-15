"""utils: Utility functions for the kmerfeatures pipeline.

author: @christinehc
"""
# imports
import numpy as np

from Bio import SeqIO


# functions
def read_fasta(fasta):
    """Read fasta file and obtain sequences, IDs, and mappings.

    Parameters
    ----------
    fasta : str
        Filename or /path/to/fasta/file.

    Returns
    -------
    (list, list)
        Tuple containing (seq_list, id_list, id2seq), where:
        seq_list : list
            List of sequences (as strings) from fasta file.
        id_list : list
            List of IDs from fasta file.

    """
    # read in sequences from the fasta file
    seq_list, id_list, = list(), list()
    with open(fasta, "r") as f:
        for record in SeqIO.parse(f, "fasta"):
            seq_list.append(str(record.seq))
            id_list.append(record.id)
    return seq_list, id_list


def read_example_index(example_index_file):
    """Read example index file and parse into dictionary.

    Parameters
    ----------
    example_index_file : str
        /path/to/example_indexfile

    Returns
    -------
    dict
        Dictionary of example indices.
        Returns empty dictionary if no index file given.

    """
    example_index = {}
    if example_index_file:
        with open(example_index_file, "r") as f:
            for line in f.readlines():
                seq_id = line.split()[0]
                example_index[seq_id] = 1.0
    return example_index


def output_to_npy(filename, id_length=6):
    """Convert feature output to numpy arrays.

    Parameters
    ----------
    filename : str
        /path/to/file.ext
    id_length : int
        Length of feature ID string (default: 6).

    Returns
    -------
    (numpy.ndarray, numpy.ndarray)
        Tuple of numpy arrays (feature_array, vector_array), where:
            feature_array : numpy.ndarray
                One-dimensional array containing feature IDs
            vector_array : np.ndarray
                N-dimensional array of vectorized feature outputs
                for all features listed in `feature_array`
        Note that the arrays are ordered; that is, the first vector
        in `vector_array` corresponds to the first ID in
        `feature_array`.

    """
    features, vectors = list(), list()
    with open(filename) as f:
        for line in f:
            line_data = line.split("\t")
            if len(line_data[0]) == id_length:
                features.append(line_data[0])
                vectors.append(np.array([float(el) for el in line_data[1:]]))
    return np.array(features), np.array(vectors)
