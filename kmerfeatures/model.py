"""model: Machine learning models for Kmer family prediction.

author: @christinehc
"""
# imports
import json
import numpy as np
import pandas as pd
from utils import get_family s

# functions
def format_data_df(filenames, label_name='family'):
    """Format Kmer sequence data into long-form dataframe.

    Dataframe consists of the following columns:
        'seq_id': Sequence IDs
        'vec': Kmer-ized vectors
        'family' (or other label name): Pre-defined label
        '{label_0, ..., label_n}': Positive identification as a
                                   given label; 0 for label, 1
                                   for non-label. Note that this
                                   gives n separate columns

    Parameters
    ----------
    filenames : type
        Description of parameter `filenames`.
    label_name : type
        Description of parameter `label_name`.

    Returns
    -------
    type
        Description of returned object.

    """
    data = {'seq_id': [], 'vec': [], label_name: []}
    for fn in filenames:
        with open(fn, 'r') as f:
            tmp = json.load(f)
            data['seq_id'] += [tmp['seq_id']]
            data['vec'] += [np.array(data['vector'])]
            data[label_name] += [get_family(fn)]

    # reformat into numpy-compatible arrays and convert to df
    data['seq_id'] = np.hstack([np.array(s) for s in data['seq_ids']])
    data['vec'] = np.concatenate([np.array(vec) for vec in data['vec']])
    data[label_name] = np.hstack(
        np.array(
            [[label] * len(vecs[i]) for i, label in enumerate(data[label_name])],
            dtype='object'
            )
        )
    data = pd.DataFrame(data)

    # get families by sizes (largest first)
    sorted_fams = data.groupby([label_name], as_index=False).count()[
        [label_name, 'seq_id']
        ].sort_values(by='seq_id', ascending=False)[label_name].values

    # define labels in binary terms (family vs. not-family)
    for topfam in sorted_fams:
        data[topfam] = [1 if fam == topfam else 0 for fam in data[label_name]]

    return data
