configfile: "config/background.yaml"

# imports
# built-in imports
import json
import pickle
from datetime import datetime
from glob import glob
from itertools import (product, repeat)
from multiprocessing import Pool
from os import makedirs
from os.path import (basename, dirname, exists, join, splitext)

# external libraries
import snekmer as skm
import numpy as np
import matplotlib.pyplot as plt
from pandas import (DataFrame, read_csv, read_json)
from Bio import SeqIO

# change matplotlib backend to non-interactive
plt.switch_backend('Agg')

# define background files
bg_files = config['input']
if isinstance(bg_files, str):
    bg_files = [bg_files]
BGS = [
    skm.utils.split_file_ext(basename(f))[0] for flist in [
        glob(join("input", bg)) for bg in bg_files
    ]
    for f in flist
]
NON_BGS, BGS = [f for f in FAS if f not in bg_files], bg_files
