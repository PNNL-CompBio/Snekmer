"""io: Input-output handling in Snekmer.

author: @christinehc

"""
# imports
import gzip
import json
import re
from os.path import basename, join

import numpy as np
import pandas as pd
from Bio import SeqIO
from .alphabet import ALPHABET_ORDER


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
        /path/to/file.ext

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
            if re.findall(r"^KMER", line_data[0]):
                kmers += line_data
    prefix = re.search(r"^KMER-[\d]+-[A-Za-z]+-", kmers[0]).group()
    return [s.strip().replace(prefix, "") for s in kmers]


def output_to_df(filename):
    """Convert kmer features and vectors into dataframe.

    Parameters
    ----------
    filename : str
        /path/to/file.txt

    Returns
    -------
    pandas.DataFrame
        DataFrame containing kmer features and vectors

    """
    features, vectors = output_to_npy(filename)
    kmers = read_output_kmers(filename)
    df = pd.DataFrame({f: v for f, v in zip(features, vectors)}).T
    df.columns = kmers
    return df.reset_index().rename(columns={"index": "seq_id"})


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
            if not re.findall(r"^KMER", line_data[0]):
                features.append(line_data[0])
                vectors.append(np.array([float(el) for el in line_data[1:]]))
    return np.array(features), np.array(vectors)


def vecfiles_to_df(files, labels=None, label_name="label"):
    """Load multiple sequence files and parse into common dataframe.

    Parameters
    ----------
    files : list
        List of (/path/to/file1.ext1, /path/to/file2.ext2, ..., etc.)
    labels : list of str or None (default: None)
        Labels for specified files; if None, excludes labels.
    label_name : str (default: 'label')
        Desired name for data column containing labels.

    Returns
    -------
    pandas.DataFrame
        Tabulated data of the form:

        | filename | seq_id | vector | vec_shape |
        |----------|--------|--------|-----------|
        | file1    | seq1   | array1 | shape1    |
        ...
        ... etc.

    """
    data = {"filename": [], "seq_id": [], "vector": [], "vec_shape": []}
    for afile in files:
        with gzip.open(afile, "rt") as f:
            # with open(afile, 'r') as f:
            tmp = json.load(f)
            data["filename"] += [basename(afile)] * len(tmp["seq_id"])
            data["seq_id"] += tmp["seq_id"]
            data["vector"] += tmp["vector"]
            data["vec_shape"] += [np.array(arr).shape for arr in tmp["vector"]]
    if str(labels) != "None":
        if len(labels) != len(data):
            raise ValueError("Number of labels must equal number of sequences.")
        data[label_name] = labels
    return pd.DataFrame(data)


def define_output_dir(alphabet, k, nested=False):
    """Create output directory name using AAR parameters.

    Parameters
    ----------
    alphabet : str or int
        Name or numerical identifier of alphabet.
        See `snekmer.alphabet.ALPHABETS` for descriptions.
    k : int
        Kmer length.
    nested : bool (default: False)
        If True, creates nested directories with alphabet, k.
        If False, returns simple output directory.

    Returns
    -------
    type
        Description of returned object.

    """
    if not nested:
        return "output"
    if not isinstance(alphabet, str):
        alphabet = ALPHABET_ORDER[alphabet]
    return join("output", alphabet, f"k-{k:02}")
