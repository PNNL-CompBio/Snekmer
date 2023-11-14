"""
Parsing methods for different file types

"""
import os, sys, csv, copy
from types import *
from SIEVEInit import *

def parse_blast_tabfile(filehandle):
    blast_results = {}

    reader = csv.reader(filehandle, delimiter="\t")
    for row in reader:
        query = row[0]
        obit = row[1].split(":")
        if len(obit) > 1:
            organism = obit[1]
        else:
            continue
        
        try:
            pid = float(row[2])/100.0
        except ValueError:
            pid = 0

        try:
            hitlen = int(row[3])
        except ValueError:
            hitlen = 0

        if not blast_results.has_key(query):
            blast_results[query] = copy.copy(ORGANISM_DICT)

        nmatches = hitlen * pid
        if nmatches > blast_results[query][organism]:
            blast_results[query][organism] = nmatches

    return blast_results
