"""io: Input-output handling for kmerfeatures pipeline.

author: @christinehc
"""
# imports
import json
import re
from os.path import basename

import numpy as np
import pandas as pd
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


def read_output_kmers(filename):
    """Extract kmer identifiers from kmerfeatures output file.

    Parameters
    ----------
    filename : str
        /path/to/output/file.txt

    Returns
    -------
    list
        Array of strings (kmer identifiers)

    """
    kmers = list()
    with open(filename) as f:
        for line in f:
            line_data = line.split("\t")
            # parse kmer outputs if detected
            if re.findall(r'^KMER', line_data[0]):
                kmers += line_data
    return [s.strip().split("-")[-1] for s in kmers]


def output_to_df(filename):
    """Short summary.

    Parameters
    ----------
    filename : type
        Description of parameter `filename`.

    Returns
    -------
    type
        Description of returned object.

    """
    features, vectors = output_to_npy(filename)
    kmers = read_output_kmers(filename)
    df = pd.DataFrame({f: v for f, v in zip(features, vectors)}).T
    df.columns = kmers
    return df.reset_index().rename(columns={'index': 'seq_id'})


def output_to_npy(filename):
    """Convert feature output to numpy arrays.

    Parameters
    ----------
    filename : str
        /path/to/output/file.txt

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
            # skip kmer outputs and only parse vectors
            if not re.findall(r'^KMER', line_data[0]):
                features.append(line_data[0])
                vectors.append(np.array([float(el) for el in line_data[1:]]))
    return np.array(features), np.array(vectors)


def vecfiles_to_df(files, labels=None, label_name='label'):
    """Load multiple sequence files and parse into common dataframe.

    Parameters
    ----------
    files : list
        Description of parameter `files`.
    labels : list or None (default: None)
        Description of parameter `labels`.
    label_name : str (default: 'label')
        Description of parameter `label_name`.

    Returns
    -------
    pandas.DataFrame
        Description of returned object.

    """
    data = {'filename': [], 'seq_id': [], 'vector': [], 'vec_shape': []}
    for afile in files:
        with open(afile, 'r') as f:
            tmp = json.load(f)
            data['filename'] += [basename(afile)] * len(tmp['seq_id'])
            data['seq_id'] += tmp['seq_id']
            data['vector'] += tmp['vector']
            data['vec_shape'] += [np.array(arr).shape for arr in tmp['vector']]
    if str(labels) != "None":
        if len(labels) != len(data):
            raise ValueError(
                'Number of labels must equal number of sequences.'
                )
        data[label_name] = labels
    return pd.DataFrame(data)
