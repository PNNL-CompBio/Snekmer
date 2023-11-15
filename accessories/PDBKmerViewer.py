#!/usr/bin/env python
"""
Standalone program to generate structures with kmer confidence score coloring.
"""
from jinja2 import Environment, FileSystemLoader
import os
import argparse

import snekmer as skm
from Bio import SeqIO

env = Environment(loader=FileSystemLoader("templates"), auto_reload=False)

def get_score_vector(fastafile, scorefile, alphabet, verbose=False):
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

    sequence_list = []
    ids_list = []
    sequence_dict = {}
    handle = open(fastafile, "rU")

    # we only need to grab the first fasta sequence
    # - actually, we should just grab it from the pdb
    #   though that can be dangerous since pdbs are notoriously
    #   bad with sequence continuity
    for record in SeqIO.parse(handle, "fasta"):
        sequence_list.append(str(record.seq))
        id = record.id
        ids_list.append(id)
        sequence_dict[id] = str(record.seq)
    handle.close()

    sequence = sequence_list[0]
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

def create_html_from_template(score_vector, pdbfile, outputfile, templatefile, verbose=None):
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
        html = template.render({"pdb_file":pdbdump,"score_vector":scoredump})
        with open(outputfile, "w") as f:
            f.write(html)

def main(fastafile=None, scorefile=None, pdbfile=None, outfile=None, alphabet=None, templatefile=None, verbose=None, **kw):
    """
    fastafile should have the first entry match with input pdbfile
    scorefile should be a tab-delimited file as output from snekmer model.
        e.g.: SSSSSVSSVSVSSV 0.23
    pdbfile should be a pdb file that matches the fastafile
    outfile the filename of the html file to output
    templatefile should be a template for the display (see examples for requirements)
    """

    score_vector = get_score_vector(fastafile, scorefile, alphabet, verbose=verbose)
    pdb_output = create_html_from_template(score_vector, pdbfile, outfile, templatefile, verbose=verbose)

OPTION_LIST = ["A program to generate features and run SIEVE models on input sequences",
                ("v", "verbose",
                 "on", None,
                 "verbose output"),
                ("f", "fastafile",
                 str, None,
                 "FASTA-format file containing protein sequences"),
                ("s", "scorefile",
                 str, None,
                 "Tab-delimited score file that contains kmer labels and weights (from RFE, e.g.)"),
                ("p", "pdbfile",
                  str, None,
                  "PDB file matching the input fasta and score"),
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
