#!/usr/bin/env python
"""
Standalone program and module to generate SIEVE feature sets from
           fasta file sequence inputs.
"""
import os, sys, gzip, random, copy, re
from types import *

from Util.Options import *
from Util.Parser import parse_blast_tabfile
from Util.SIEVEInit import *

from distutils.spawn import find_executable

#try:
from Bio import SeqIO
#from Bio.Blast import NCBIStandalone
#from Bio.Blast import NCBIXML

#from Bio.Clustalw import MultipleAlignCL
#from Bio import Clustalw
from Bio.Align import AlignInfo
from Bio.Alphabet import IUPAC

#except ImportError:
#    sys.stderr.write("BioPython not installed correctly (see http://biopython.org)\n")
#    sys.exit(-1)

OPTION_LIST = ["A program to generate features and run SIEVE models on input sequences",
                "None",
                (None, "v", "verbose",
                 "on", None, None, None,
                 "verbose output"),
                (None, "f", "fastafile",
                 "str", None, None, None,
                 "FASTA-format file containing protein sequences"),
                (None, "i", "indexfile",
                 "str", None, None, None,
                 "text file with the first column containing ids of positive examples"),
                (None, "s", "scorefile",
                 "str", None, None, None,
                 "Tab-delimited score file that contains kmer labels and weights (from RFE, e.g.)"),
                (None, "g", "gap",
                 "int", None, None, None,
                 "Gap length- assumed to be center of the kmer"),
                (None, "G", "giids",
                 "int", None, None, None,
                 "Parse fasta id line as gi line (split on '|')"),
                (None, "o", "outfile",
                 "str", None, None, None,
                 "Output HTML file with annotated sequences")
    ]

def baseconvert(n, k=None, digits=None, **kw):
    """convert positive decimal integer n to equivalent in another """

    #digits = "0123456789abcdefghijklmnopqrstuvwxyz"
    digits = digits or "ACDEFGHIKLMNPQRSTVWY"
    base = len(digits)

    try:
        n = int(n)
        base = int(base)
    except:
        return ""

    if n < 0 or base < 2 or base > 36:
        return ""

    s = ""
    while 1:
        r = n % base
        s = digits[r] + s
        n = n / base
        if n == 0:
            break

    if len(s) < k:
        for i in range(len(s), k):
            s = "%s%s" % (digits[0], s)

    return s

def main(fastafile=None, scorefile=None, outfile=None, indexfile=None, giids=True, gap=None, verbose=None, **kw):
    """
    fastafile should have sequences that you want to 'annotate' by highlighting important kmers
    scorefile should be a tab-delimited file as output from RFE.
        e.g.: KMER-2-NAT-CG   7       19.1355 1
        kmer label [tab] rfe iteration [tab] score [tab] keep/discard (unused here)
    outfile the filename of the html file to output

    """
    fd = open(scorefile, "r")
    scores = fd.readlines()
    fd.close()
    score_dict = {}
    reg_dict = {}

    for line in scores:
        bits = line.split("\t")
        label = bits[0]
        score = float(bits[2])
        kmer = label.split("-")[-1]
        alpha = label.split("-")[-2]
        score_dict[kmer] = score

        alpha = alpha or "NAT"
        if alpha != "NAT":
            alphadict = ALPHABETCODE[alpha]

            reg = ""
            for bit in kmer:
                reg = "%s[%s]" % (reg, alphadict[bit])
            reg_dict[kmer] = reg

    posids = []
    if indexfile:
        fd = open(indexfile, "r")
        ids = fd.readlines()
        fd.close()

        for line in ids:
            bits = line.split("\t")
            label = bits[0]
            posids.append(label)

    sequence_list = []
    ids_list = []
    sequence_dict = {}
    handle = open(fastafile, "rU")
    for record in SeqIO.parse(handle, "fasta"):
        sequence_list.append(str(record.seq))
        id = record.id
        if giids:
            id = id.split("|")[1]
        ids_list.append(id)
        sequence_dict[id] = str(record.seq)
    handle.close()

    if indexfile:
        # reorder the ids_list
        pos_ids = []
        neg_ids = []
        for id in ids_list:
            if id in posids:
                pos_ids.append(id)
            else:
                neg_ids.append(id)
        ids_list = pos_ids + neg_ids

    out = "<pre><!-- body {font-family:monospace}--><table>"
    for id in ids_list:
        this_line = "<tr><td>%s</td><td>" % id
        if id in posids:
            this_line = "<tr><td><font color='red'>%s</font></td><td>" % id

        seq = sequence_dict[id]

        html_seq = ""

        if gap:
            gap_dict = {}
            for i in range(1, len(seq)):
                for k in score_dict.keys():
                    l = len(k)/2
                    ke = k[l:]
                    ks = k[:l]
                    #print k, ks, ke, seq[i:i+l], seq[i+l+gap:i+l+l+gap], i, i+l, i+l+gap, i+l+l+gap
                    if seq[i:(i+l)] == ks and seq[(i+gap):(i+l+gap)] == ke:
                        print("matched!", k, i, i+gap)
                        score = score_dict[k]
                        for j in range(0,l):
                            gap_dict[i+j] = score > 0
                            gap_dict[i+j+gap] = score > 0
            for i in range(1, len(seq)):
                c = seq[i]
                if gap_dict.has_key(i):
                    color = "green"
                    if gap_dict[i]:
                        color = "red"
                    c = "<font color=%s>%s</font>" % (color, c)
                html_seq += c
        else:
            i = 1
            done = False
            while not done:
                score = 0

                for k in score_dict.keys():
                    l = len(k)

                    if reg_dict and re.match(reg_dict[k], seq[i:i+l]):
                        score = score_dict[k]
                        s = seq[i:i+l]
                        break

                    elif seq[i:i+l] == k:
                        score = score_dict[k]
                        s = seq[i:i+l]
                        break

                if score != 0:
                    if score > 0:
                        color = "red"
                    else:
                        color = "green"
                    i += l
                    print("%s\t%s\t%s\t%s\t%d" % (s, score, id, k, i))
                    html_seq = "%s<font color=%s>%s</font>" % (html_seq, color, s)
                else:
                    if i < len(seq):
                        html_seq = "%s%s" % (html_seq, seq[i])

                    i += 1

                if i > len(seq):
                    done = True

        this_line = "%s%s" % (this_line, html_seq)
        this_line = "%s</td></tr>\n" % this_line
        out = "%s%s" % (out, this_line)

    out = "%s</table></font></pre>" % out
    #print out

    fd = open(outfile, "w")
    fd.write(out)
    fd.close()

if __name__ == "__main__":
    optdict, infiles = process_options(OPTION_LIST)
    apply(main, infiles, optdict)
