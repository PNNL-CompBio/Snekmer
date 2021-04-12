"""utils: Utility functions for the kmerfeatures pipeline.

author: @christinehc
"""
# imports
from datetime import datetime
import re
from os.path import basename

import numpy as np
import pandas as pd
from Bio import SeqIO


# functions
def count_sequences(fasta):
    """Count the number of sequences contained in a fasta file.

    Parameters
    ----------
    fasta : str
        Filename or /path/to/fasta/file.

    Returns
    -------
    int
        Number of sequences in the input file.

    """
    return len([record for record in SeqIO.parse(fasta, "fasta")])


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

    Note: If no family matching the given regular expression is
    found, the original filename is returned, stripped of any file
    extension (i.e. the file's basename).

    Parameters
    ----------
    filename : str
        path/to/filename.ext
    regex : str or r-string
        Regular expression for matching a family name
        (default: "[a-z]{3}[A-Z]{1}").
        The default expression is 3 lowercase letters followed
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
    filename = basename(filename)
    s = '_'.join(
        filename.split('.')[:-1]
        ).replace('-', '_').replace(' ', '_')
    # modify regex to only search between underscores
    regex = r'(?<=_)' + r'{}'.format(regex) + r'(?=_)'

    # return list output
    if not return_first:
        re.findall(regex, s)

    # return string output
    if re.search(regex, s):
        return re.search(regex, s).group()
    return filename.split('.')[0]


def format_timedelta(timedelta):
    """Format datetime.timedelta object as str with units.

    Parameters
    ----------
    timedelta : datetime.timedelta object
        A datetime.timedelta object.

    Returns
    -------
    str
        Formatted string with the following format:
            %%h %%m %%.%%s
            (e.g. '01h 41m 20.01s')

    """
    # get milliseconds and round
    milliseconds = round(timedelta.microseconds / 1000)

    # get H, M, S
    minutes, seconds = divmod(timedelta.seconds, 60)
    hours, minutes = divmod(minutes, 60)

    return f"{hours}h {minutes:02}m {seconds:02}.{milliseconds}s"
