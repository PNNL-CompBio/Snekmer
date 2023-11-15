#!/usr/bin/env python
"""
Standalone program to generate structures with kmer confidence score coloring.
"""
from jinja2 import Environment, FileSystemLoader
import os
import argparse

import snekmer as skm
from Bio import SeqIO
from Bio.PDB import *

env = Environment(loader=FileSystemLoader("templates"), auto_reload=False)

def get_score_vector(scorefile, alphabet, sequence=None, verbose=False):
    fd = open(scorefile, "r")
    scores = fd.readlines()
    fd.close()
    score_dict = {}
    reg_dict = {}

    for line in scores:
        bits = line.split(",")
        if not bits[0] == 'kmer':
            kmer = bits[0]
            score = float(bits[1])
            #kmer = label.split("-")[-1]
            #alpha = label.split("-")[-2]
            score_dict[kmer] = score

    # this may be unwise if we wanted to do something creative
    #  with varibale k lengths
    k = len(kmer)

    sequence = skm.vectorize.reduce(sequence, alphabet)

    # create an emptyp scoring vector for the sequence
    seqscore = [-1.0,]*len(sequence)

    for i in range(0,len(sequence)):
        seqkey = sequence[i:i+k]
        if seqkey in score_dict:
            score = score_dict[seqkey]
            for j in range(0,k):
                if seqscore[i+j] < score:
                    # maximum score for each position
                    seqscore[i+j] = score

    if verbose:
        for i in range(0,len(sequence)):
            print("%d\t%s\t%.3f" % (i, sequence[i], seqscore[i]))

    return(seqscore)

def get_sequence_from_pdb(pdbfile):
    d3to1 = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
            'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
            'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
            'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

    # run parser
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('struct', pdbfile)

    # iterate each model, chain, and residue
    # printing out the sequence for each chain
    seqs = {}
    for model in structure:
        for chain in model:
            id = "%s:%s" % (model.id, chain.id)
            seq = []
            pdbstart = -1
            for residue in chain:
                if pdbstart < 0:
                    pdbstart = residue.get_id()[1]
                # warning: ignores unknown codes
                if residue.resname in d3to1:
                    seq.append(d3to1[residue.resname])
            seqs[id] = (pdbstart, "".join(seq))

    return(seqs)

def create_html_from_template(score_vector, pdbfile, outputfile,
                                templatefile, pdbstart=1, verbose=None):
        """Create an HTML 3Dmol viewer for the pdb and kmer scores.

        Parameters
        ----------
        score_vector : list
            List of scores for each residue in the pdb file
        pdbfile : str
            Filename of the input PDB file
        templatefile : str
            templatefile to use as a base
        outputfile : str
            Filename of the output html 3Dmol viewer

        Returns
        -------
        None
            Creates file and exits.

        """
        # just slurp it all in at once
        pdbdump = open(pdbfile, "r").read()
        scoredump = str(score_vector)

        template = env.get_template(templatefile)
        html = template.render({"pdb_file":pdbdump,
                                "score_vector":scoredump,
                                "pdb_start":pdbstart})
        with open(outputfile, "w") as f:
            f.write(html)

def get_pdbfile_from_id(pdbid):
    pdbl = PDBList()
    pdbfile = pdbl.retrieve_pdb_file(pdbid, file_format='pdb', pdir='pdb', obsolete=False)

    return(pdbfile)

def main(scorefile=None, pdbfile=None, pdbid=None, outfile=None, alphabet=None, templatefile=None, verbose=None, **kw):
    """
    scorefile should be a tab-delimited file as output from snekmer model.
        e.g.: SSSSSVSSVSVSSV 0.23
    pdbfile should be a pdb file to view
    outfile the filename of the html file to output
    templatefile should be a template for the display (see examples for requirements)
    """

    # get pdbfile from the PDB and save it locally
    # FIXME: this dumps the file in the current directory
    if pdbid:
        pdbfile = get_pdbfile_from_id(pdbid)

    sequences = get_sequence_from_pdb(pdbfile)

    # this will only grab the first record from the pdb and
    #     not sure this id will be applicable to all pdbs
    (pdbstart, sequence) = sequences["0:A"]

    score_vector = get_score_vector(scorefile, alphabet, sequence=sequence, verbose=verbose)
    pdb_output = create_html_from_template(score_vector, pdbfile, outfile,
                                            templatefile, pdbstart=pdbstart,
                                            verbose=verbose)

OPTION_LIST = ["A program to generate features and run SIEVE models on input sequences",
                ("v", "verbose",
                 "on", None,
                 "verbose output"),
                ("s", "scorefile",
                 str, None,
                 "Tab-delimited score file that contains kmer labels and weights (from RFE, e.g.)"),
                ("p", "pdbfile",
                  str, None,
                  "PDB file to view"),
                ("P", "pdbid",
                  str, None,
                  "Valid identifier of a PDB record to retrieve"),
                ("o", "outfile",
                 str, "output.html",
                 "Output HTML file with annotated sequences"),
                 ("t", "templatefile",
                 str, "pdb_package_template.html",
                 "HTML template file for 3Dmol viewer"),
                 ("a", "alphabet",
                 str, "NAT",
                 "Alphabet for kmer encoding")
    ]

def get_argument_parser():
    parser = {}

    parser["main"] = argparse.ArgumentParser(
        description=(OPTION_LIST[0])
    )

    for optchunk in OPTION_LIST[1:]:
        thistype = optchunk[2]
        thesekw = {"type":thistype}
        action = 'store'
        if thistype == "on":
            thesekw = {}
            action = 'store_true'

        parser["main"].add_argument(
            ("--%s" % optchunk[1]),
            ("-%s" % optchunk[0]),
            dest=optchunk[1],
            action=action,
            default=optchunk[3],
            help=optchunk[4],
            **thesekw
            )
    return parser

if __name__ == "__main__":
    parser = get_argument_parser()

    # parse args
    args = parser["main"].parse_args()

    main(**vars(args))
