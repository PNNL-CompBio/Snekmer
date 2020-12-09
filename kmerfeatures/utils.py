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


def read_example_index(example_indexfile):
    """Read example index file and parse into dictionary.

    Parameters
    ----------
    example_indexfile : str
        /path/to/example_indexfile

    Returns
    -------
    dict
        Dictionary of example indices.
        Returns empty dictionary if no index file given.

    """
    example_index = {}
    if example_indexfile:
        with open(example_indexfile, "r") as f:
            for line in f.readlines():
                seq_id = line.split()[0]
                example_index[seq_id] = 1.0
    return example_index
