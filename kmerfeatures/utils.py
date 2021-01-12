"""utils: Utility functions for the kmerfeatures pipeline.

author: @christinehc
"""
# imports
import re

import numpy as np
import pandas as pd

from Bio import SeqIO
from os.path import basename


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


def get_output_kmers(filename):
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
    kmers = get_output_kmers(filename)
    df = pd.DataFrame({f: v for f, v in zip(features, vectors)}).T
    df.columns = kmers
    return df.reset_index().rename(columns={'index': 'seq_id'})


def parse_fasta_description(fasta, df=True):
    """Parse description field of a FASTA file.

    Parameters
    ----------
    fasta : str
        path/to/fasta.file
    df : bool
        If True, returns output as a pandas DataFrame (default: True)
        If False, returns output as a dictionary.

    Returns
    -------
    dict or pandas.DataFrame
        Parsed values as {key: value} pairs or collated in DataFrame.

    """
    # collect parsed description data into dict
    pdict = dict()
    with open(fasta, 'r') as f:
        for record in SeqIO.parse(f, 'fasta'):
            pdict[record.id] = dict()
            pdict[record.id]['protein_family'] = get_family(fasta)
            s = f"{record.description} "  # trailing space needed for regex
            parsed = re.findall(r'([\w]+[=][\w\", ]+)(?= )(?!=)', s)
            for item in parsed:
                key = item.split("=")[0].lower()
                val = item.split("=")[1].replace('"', '')
    #             print(key, val)
                i = 1
                while key in pdict[record.id].keys():
                    key = f"{key.rstrip('0123456789_')}_{i}"
                    i += 1
                pdict[record.id][key] = val
            pdict[record.id]['filename'] = basename(fasta)
    if df:
        pdict = pd.DataFrame(pdict).T.reset_index().rename(
            columns={'index': 'sequence_id'}
            )
    return pdict


def get_family(filename, regex='[a-z]{3}[A-Z]{1}', return_first=True):
    """Extract family from filename given regular expression format.

    Parameters
    ----------
    filename : str
        path/to/filename.ext
    regex : str or r-string
        Regular expression for matching a family name
        (default: "[a-z]{3}[A-Z]{1}").
        The default expression looks for 3 lowercase letters followed
            by one uppercase letter. To write a custom regular
            expression, see https://docs.python.org/3/library/re.html
            for more details on using the built-in re library.
    return_first : bool
        If True, returns only the first occurrence (default: True).
        If False, returns a list of all occurring regex matches.

    Returns
    -------
    str or list
        Family name or names

    """
    # extract and simplify file basename
    s = '_'.join(
        basename(filename).split('.')[:-1]
        ).replace('-', '_').replace(' ', '_')
    # modify regex to only search between underscores
    regex = r'(?<=_)' + r'{}'.format(regex) + r'(?=_)'
    if return_first:
        return re.search(regex, s).group()
    return re.findall(regex, s)
