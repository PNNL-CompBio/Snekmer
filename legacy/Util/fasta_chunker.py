#!/usr/bin/env python
"""
Standalone program to run SIEVE on different models from user input sequences.
"""
import os, sys, gzip
from types import *

from Options import *
from Bio import SeqIO

OPTION_LIST = ["A program to break fasta files into chunks",
                None,
                (None, "v", "verbose",
                 "on", None, None, None,
                 "verbose output"),
                (None, "f", "fastafile",
                 "str", None, None, None,
                 "FASTA-format file containing protein sequences"),
                (None, "c", "chunk_size",
                 "int", None, None, None,
                 "Number of proteins to keep in each chunk"),
                (None, "o", "output_filebase",
                 "str", None, None, None,
                 "Base filename for output chunks")]

def main(fastafile=None, chunk_size=None, output_filebase=None, verbose=None, **kw):
    if output_filebase == None:
        output_filebase = os.path.splitext(fastafile)[0]

    chunk_size = chunk_size or 1000
    
    handle = open(fastafile, "rU")
    chunk = []
    nprots = 0
    chunk_n = 0
    for prot in SeqIO.parse(handle, "fasta"):
        if not divmod(nprots, chunk_size)[1] and chunk:
            chunk_n += 1
            fname = "%s_%04d.fasta" % (output_filebase, chunk_n)
            cfile = open(fname, "w")
            for c in chunk:
                cfile.write(">%s %s\n%s\n" % (c.id, c.description, str(c.seq)))
            cfile.close()
            chunk = []
        chunk.append(prot)
        nprots += 1
        if verbose:
            print "Protein %d going to chunk %d" % (nprots, chunk_n)

    chunk_n += 1
    fname = "%s_%04d.fasta" % (output_filebase, chunk_n)
    cfile = open(fname, "w")
    for c in chunk:
        cfile.write(">%s %s\n%s\n" % (c.id, c.description, str(c.seq)))
    cfile.close()
            
    handle.close()

if __name__ == "__main__":
    optdict, infiles = process_options(OPTION_LIST)
    apply(main, infiles, optdict)
