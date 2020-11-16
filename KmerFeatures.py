#!/usr/bin/env python
"""
Standalone program and module to generate SIEVE feature sets from
           fasta file sequence inputs.
"""
import os, sys, gzip, random, copy
from types import *

from Util.Options import *
from Util.SIEVEInit import *

try:
    from Bio import SeqIO
    #from Bio.Alphabet import IUPAC
    
except ImportError:
    sys.stderr.write("BioPython not installed correctly (see http://biopython.org)\n")
    sys.exit(-1)

OPTION_LIST = ["A program to generate features and run SIEVE models on input sequences",
                "None",
                (None, "v", "verbose",
                 "on", None, None, None,
                 "verbose output"),
                (None, "p", "parameterfile",
                 "str", None, None, None,
                 "Filename of a parameter file to use"),
                (None, "f", "fastafile",
                 "str", None, None, None,
                 "FASTA-format file containing protein sequences"),
                (None, "k", "kmer",
                 "int", 3, None, None,
                 "kmer length for features"),
                (None, "s", "start",
                  "int", None, None, None,
                 "first residue of sequences to start generating features"),
                (None, "n", "end",
                 "int", None, None, None,
                 "last residue of sequences to start generating features"),
                (None, "D", "nucleotide",
                 "off", None, None, None,
                 "process nucleotide sequences"),
                (None, "M", "map_function",
                 "str", None, ("reduced_alphabet_0", "reduced_alphabet_1", "reduced_alphabet_2", "reduced_alphabet_3", "reduced_alphabet_4"),
                 ("Hydrophobic/hydrophilic", "Standard 7", "Solvent accessibility", "LoveHateCharge", "LoveHateBadstruct"),
                 "mapping function for reduced amino acid alphabets"),
                (None, "r", "randomize_alphabet",
                 "off", None, None, None,
                 "Randomize alphabet used - but using the same number of categories and distribution as specificed by the map_function"),
                (None, "R", "min_rep_thresh",
                 "float", 1, None, None,
                 "minimum number of sequences to include feature for prefiltering. 0>R<1 is percentage of input sequences"),
                (None, "e", "example_indexfile",
                 "str", None, None, None,
                 "File containing identifiers of positive examples for sieve output records or gist .class files"),
                (None, "m", "features_output_format",
                 "str", "simple", ("simple", "gist", "sieve"), ("Simple tab-delimited format", "Input for gist", "Input for sieve"),
                 "Format for outputting feature sets"),
                (None, "o", "features_output_filebase",
                 "str", None, None, None,
                 "Filename base (no suffix) for output features"),
                (None, "D", "filter_duplicates",
                 "on", None, None, None,
                 "Filter out duplicate identifiers from the fasta files and/or blast output"),
                (None, "Q", "output_shuffled_sequences",
                 "int", None, None, None,
                 "Output shuffled sequences in a fasta file for each input sequence"),
                (None, "F", "feature_set",
                 "str", None, None, None,
                 "Specify which features to include in the model"),
                (None, "K", "kmer_output",
                "on", None, None, None,
                "Verbose output of kmer positions for each sequence [under development]"),
                (None, "w", "walk",
                 "off", None, None, None,
                 "Perform a kmer walk on the input fasta file to get an idea of the kmer representation")
    ]


def output_features(feature_sets=None, format=None, output_filename=None, labels=None, write_mode="w", **kw):
    if format == "gist" or format == "both":
        trainout = "%s.train" % output_filename
        classout = "%s.class" % output_filename
        trainhandle = open(trainout, write_mode)

        trainhandle.write("corner")
        if labels:
            for label in labels:
                trainhandle.write("\t%s" % label)
        else:
            for i in range(len(feature_sets[0])-1):
                trainhandle.write("\tlabel%d" % i)
        trainhandle.write("\n")
        
        classhandle = open(classout, write_mode)
        classhandle.write("corner\tclass\n")
        
        for features in feature_sets:
            output_gist_features(features=features, filehandle=trainhandle, **kw)
            output_gist_class(features=features, filehandle=classhandle, **kw)
        trainhandle.close()
        classhandle.close()
            
    if format == "sieve" or format == "both":
        patternout = "%s.pattern" % output_filename
        outhandle = open(patternout, write_mode)
        for features in feature_sets:
            output_sieve_features(features=features, filehandle=outhandle, **kw)
        outhandle.close()

    if format == "matrix":
        fileout = "%s.txt" % output_filename
        filehandle = open(fileout, write_mode)
        if labels:
            filehandle.write("%s" % labels[0])
            for label in labels[1:]:
                filehandle.write("\t%s" % label)
            filehandle.write("\n")
        if feature_sets:
            for features in feature_sets:
                output_gist_features(features=features, filehandle=filehandle, **kw)
        filehandle.close()

def output_gist_features(features=None, filehandle=None, **kw):
    filehandle.write("%s" % features[0])
    for f in features[1:]:
        filehandle.write("\t%s" % f)
    filehandle.write("\n")
    filehandle.flush()

def output_gist_class(features=None, filehandle=None, example_index={}, **kw):
    id = features[0]
    value = example_index.get(id, -1)
    filehandle.write("%s\t%d\n" % (features[0], value))
    filehandle.flush()

def output_sieve_features(features=None, filehandle=None, example_index={}, **kw):
    id = features[0]
    value = example_index.get(id, 0.0)
    filehandle.write("pattern\t%s\t%d\n" % (id, len(features)-1))
    filehandle.write("\tinput\t%s" % id)
    for f in features[1:]:
        filehandle.write("\t%s" % f)
    filehandle.write("\n")
    filehandle.write("\toutput\t%s\t%d\n" % (id, value))
    filehandle.flush()
    
def string_vectorize(residues=None, sequence=None, kmer=3, start=None, end=None, map_function=None, return_labels=None,
                     feature_dict=None, filter_list=None, exclusion_list=None, return_dict=None, kmer_output=None, **kw):
    def identity(character, **kw):
        return character

    def reduce_alphabet(character, mapping=None, **kw):
        for key in mapping.keys():
            if character in key:
                return mapping[key]

    residues = residues or StandardAlphabet
    map_function = map_function or identity
    mapping = None
    map_name = "NAT"

    # this is specifically for when we create a random alphabet
    #   that we want to apply to a lot of sequences
    if type(map_function) == list:
        
        kw["mapping"] = map_function[2]
        residues = map_function[0]
        map_name = map_function[1]
        map_function = reduce_alphabet

    if map_function == "reduced_alphabet_0":
        map_function = reduce_alphabet
        kw["mapping"] = get_alphabets()["reduced_alphabet_0"]
        residues = "SV"
        map_name = "RED0"
        #residues = reduced_alphabet_0.values()

    if map_function == "reduced_alphabet_1":
        map_function = reduce_alphabet
        kw["mapping"] = get_alphabets()["reduced_alphabet_1"]
        residues = "APFNDKC"
        map_name = "RED1"
        #residues = reduced_alphabet_1.values()

    if map_function == "reduced_alphabet_2":
        map_function = reduce_alphabet
        kw["mapping"] = get_alphabets()["reduced_alphabet_2"]
        residues = "CAP"
        map_name = "RED2"
        #residues = reduced_alphabet_2.values()

    if map_function == "reduced_alphabet_3":
        map_function = reduce_alphabet
        kw["mapping"] = get_alphabets()["reduced_alphabet_3"]
        residues = "LHC"
        map_name = "RED3"

    if map_function == "reduced_alphabet_4":
        map_function = reduce_alphabet
        kw["mapping"] = get_alphabets()["reduced_alphabet_4"]
        residues = "LHB"
        map_name = "RED4"

    # we should provide the ability to include a bunch of strings that can be used to vectorize
    #if pow(len(residues), kmer) > 20000:
    #    sys.stderr.write("Error: this will generate more than 20k inputs\n")
    #    return

    if return_labels:
        labels = []
        if filter_list:
            for this in filter_list:
                label = "KMER-%d-%s-%s" % (kmer, map_name, this)
                labels.append(label)
            return labels

        else:
            for bit in range(pow(len(residues), kmer)):
                label = "KMER-%d-%s-%s" % (kmer, map_name, baseconvert(bit, k=kmer, digits=residues))
                labels.append(label)
            return labels

    tstart = start or 0
    if tstart < 0:
        tstart = len(sequence) + tstart
        if tstart < 0:
            tstart = 0

    tend = end or len(sequence)-kmer

    results = {}
    if feature_dict:
        results = feature_dict

    elif filter_list:

        # there is a more elegant way of doing this- and probably clever too
        for item in filter_list:
            results[item] = 0

    for i in range(tstart, tend):
        kmap = sequence[i:i+kmer]

        # do mapping to a reduced alphabet, e.g.
        kstring = ""
        for c in kmap:
            cc = map_function(c, **kw)
            
            # this has the effect of omitting unrecognized characters, which may be undesireable in some cases
            if cc == None:
                continue
            kstring += cc

        if len(kstring) < kmer:
            # this happens when there are unrecognized characters
            # and we need to not include these
            continue

        # if we don't find the kstring in the filter_list
        #   then we skip.
        #if filter_dict:
        #    print(filter_dict.keys())
        #    print(kstring)
        
        if kmer_output:
          print(i, "\t", kmap, "\t", kstring, "\t1", )
        
        if filter_list and not kstring in filter_list:
            #print("hoopla")
            continue
        
        #FILTER HERE
        
        # a list of strings to exclude from consideration
        #   this is from the kmer_walk approach and should
        #   be shorter sequences (though could also be the
        #   same length)
        if exclusion_list:
            breaker = 0
            for bit in exclusion_list:
                if kstring.find(bit) > 0:
                    breaker = 1
                    break
            if breaker: next
            
        if not kstring in results:
            results[kstring] = 0
        results[kstring] += 1

    if return_dict:
        return results

    if filter_list:
        results_ordered = []
        for item in filter_list:
            results_ordered.append(results[item])
        return results_ordered
    
    return results.values()

    # old stuff below
    results_vector = []
    for j in range(0, pow(len(residues), kmer)):
        results_vector.append(0.0)

    for key in results.keys():
        # base conversion- though we're using base 20 (A-Y, minus B, J, O, U) and not base 25 (A-Y)
        # for the standard aa alphabet. We can also use various reduced alphabets
        x = 0
        
        if len(key) < kmer:
            continue

        for k in range(kmer):
            # provide a unique index for all possible strings
            x += residues.find(key[k]) * pow(len(residues), k)
        results_vector[x] = results[key]

    return results_vector

def scramble_sequence(id=None, sequence=None, n=None, no_id_modifier=None, first_residue_special=True, example_index=None, **kw):
    # takes a sequence and an identifier and returns lists of n scrambled sequences
    #       and numbered identifiers
    seqlist = []
    idlist = []

    seq = []
    # protect the N-terminal M, e.g.
    if first_residue_special:
        start_residue = sequence[0]
        sequence = sequence[1:]
    else:
        start_residue = ""
    
    for c in sequence:
        seq.append(c)

    for i in range(n):
        random.shuffle(seq)
        thisseq = start_residue
        for c in seq:
            thisseq += c
            
        seqlist.append(thisseq)
        if no_id_modifier:
            idlist.append(id)
        else:
            sid = "%s_shuffle_%d" % (id, i)
            idlist.append(sid)
        example_index[sid] = -1.0
    return idlist, seqlist, example_index

def make_n_terminal_fusions(id=None, filename=None, **kw):
    sequence_list = []
    id_list = []
    handle = open(filename, "r")
    for record in SeqIO.parse(handle, "fasta"):
        thisid = "%s-%s" % (record.id, id)
        id_list.append(thisid)
        sequence_list.append(str(record.seq))
    handle.close()

    return id_list, sequence_list

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
        n = int(n)
        if n == 0:
            break

    if len(s) < k:
        for i in range(len(s), k):
            s = "%s%s" % (digits[0], s)

    return s

def kmer_walk(fastafile=None, maxk=20, seq_dump=False, **kw):
    # next read in sequences from the fasta file
    sequence_list = []
    sequence_dict = {}
    ids_list = []
    handle = open(fastafile, "r")
    for record in SeqIO.parse(handle, "fasta"):
        sequence_list.append(str(record.seq))
        ids_list.append(record.id)            
        sequence_dict[record.id] = str(record.seq)
    handle.close()

    # need to make this a percentage - eliminate the bottom XX% of kmers by repreesentation
    minthresh = 10

    exclusion_list = None
    for kmer in range(1, maxk):
        # build a feature dict for this kmer and these sequences
        feature_dict = {}
        
        # short circuit the exclusion list (very slow) but allow it to count
        exclusion_list = []
        for i in range(len(sequence_list)):
            sequence = sequence_list[i]
            id = ids_list[i]

            feature_dict = string_vectorize(sequence=sequence, kmer=kmer, map_function="reduced_alphabet_0", feature_dict=feature_dict,
                                             exclusion_list=exclusion_list, return_dict=True)

        exclusion_list = []
        for key in feature_dict.keys():
            if feature_dict[key] < minthresh:
                exclusion_list.append(key)
        print("Kmer %d, number of remaining features %d total, number of remaining features occuring more than %d times %d, of %g possible, %g%%" %
              (kmer, len(feature_dict.keys()), minthresh, len(feature_dict.keys())-len(exclusion_list), 20**kmer, (len(feature_dict.keys())-len(exclusion_list))/20**kmer))

    if seq_dump:
        for seq in feature_dict.keys():
            if feature_dict[seq] > minthresh:
                print(seq)

def define_feature_space(sequence_dict=None, kmer=None, map_function=None, start=None, end=None, residues=None, min_rep_thresh=2, **kw):
    # this routine will return

    feature_dict = {}
    
    for id, seq in sequence_dict.items():
        feature_dict = string_vectorize(sequence=seq, kmer=kmer, map_function=map_function, feature_dict=feature_dict,
                                         start=start, end=end, residues=residues, return_dict=True)

    # if this is between 0 and 1 then it's a percentage
    if min_rep_thresh < 1 and min_rep_thresh > 0:
        min_rep_thresh = len(feature_dict.keys()) * min_rep_thresh
        
    # filter out all those below the min_rep_thresh
    if min_rep_thresh:
        filter_dict = {}
        for key in feature_dict.keys():
            if feature_dict[key] >= min_rep_thresh:
                filter_dict[key] = feature_dict[key]
    else:
        filter_dict = feature_dict
    
    return filter_dict
        

def main(fastafile=None, example_indexfile=None, features_output_format=None, features_output_filebase=None,
         filter_duplicates=None, shuffle_n=None, output_shuffled_sequences=None, n_terminal_file=None,
         svm_ism=None, feature_set=None, kmer=None, start=None, end=None, nucleotide=None, kmer_output=None, 
         map_function=None, min_rep_thresh=None, randomize_alphabet=0, walk=None, verbose=None, **kw):

    #randomize_alphabet = 1
    
    if walk:
        return kmer_walk(fastafile=fastafile)

    # next read in sequences from the fasta file
    sequence_list = []
    sequence_dict = {}
    ids_list = []
    handle = open(fastafile, "r")
    for record in SeqIO.parse(handle, "fasta"):
        sequence_list.append(str(record.seq))
        ids_list.append(record.id)            
        sequence_dict[record.id] = str(record.seq)
    handle.close()
    
    # this indexfile (optional) contains the identifiers of positive examples for feature output
    example_index = {}
    if example_indexfile:
        eind = open(example_indexfile, "r")
        for line in eind.readlines():
            id = line.split()[0]
            example_index[id] = 1.0

    feature_sets = []

    if randomize_alphabet:
        this = []
        
        alpha = get_alphabets()[map_function]
        
        residues = alpha["_keys"]
        map_name = "RND%s" % map_function[-1]

        rand_alphabet = {}

        # this gives us a string that doesn't repeat to use as the new keys
        randstr = ''.join(random.sample("ACDEFGHIKLMNPQRSTVWY", 20))
        for keyi in alpha.keys():
            if keyi == "_key":
                rand_alphabet["_key"] = alpha["_key"]
                continue
            key = randstr[0:len(keyi)]

            # trim off this sequence
            randstr = randstr[len(keyi):]
            
            rand_alphabet[key] = alpha[keyi]

        map_function = [residues, map_name, rand_alphabet]

    if features_output_filebase == None:
        features_output_filebase = fastafile

    if feature_set == None:
        # NEW: prefilter the entire fasta 
        # this may take awhile but will cut down on the size of the feature set by a lot
        filter_dict = define_feature_space(sequence_dict=sequence_dict, kmer=kmer,
                                        map_function=map_function, start=start, end=end,
                                        min_rep_thresh=min_rep_thresh)

        if verbose:
            print("Feature space: %d kmers with more than %d representation in %d sequences" % (len(filter_dict.keys()), min_rep_thresh, len(sequence_dict.keys())))
            
        #print(filter_dict.keys())
        filter_list = filter_dict.keys()
        if len(filter_list) == 0:
            print("Terminating")
            return

    else:
        # we read in a list of feature ids to use from a file
        # NOTE: not doing any check on the format of these ids
        filter_list = []
        handle = open(feature_set, "r")
        for line in handle.readlines():
              filter_list.append(line.split()[0])
        handle.close()
    
    # this will catch multiples of the same identifier- but may result in problems if there's more than
    #      one of the same identifier (i.e. user's files might be messy this way)
    seen = []

    i = 0
    first = True
    for i in range(len(sequence_list)):
        sequence = sequence_list[i]
        id = ids_list[i]

        if filter_duplicates and id in seen:
            continue
        seen.append(id)

        sequences = [sequence,]
        ids = [id,]


        # shuffle the N-terminal sequence N times and run SIEVE on the wt and each shuffled sequence
        if shuffle_n:
            example_index[id] = 1.0
            scid_list, scramble_list, example_index = scramble_sequence(id=id, sequence=sequence[:30], n=shuffle_n, example_index=example_index)
            sequences += scramble_list
            ids += scid_list

            if output_shuffled_sequences:
                filename = "%s_shuffled.fasta" % id
                fh = open(filename, "w")
                for i in range(len(ids)):
                    id = ids[i]
                    seq = sequences[i]
                    fh.write(">%s\n%s\n" % (id, seq))
                fh.close()

        if n_terminal_file:
            addid_list, addseq_list = make_n_terminal_fusions(id=id, filename=n_terminal_file)
            sequences += addseq_list
            ids += addid_list

        residues = None
        if nucleotide:
            residues = "ACGT"


        labels = []
        for j in range(len(sequences)):
            sequence = sequences[j]
            id = ids[j]

            if verbose:
                print("Constructing features for sequence %s" % id)

            features = [id,]
            
            features += string_vectorize(sequence=sequence, kmer=kmer, start=start, end=end, map_function=map_function, residues=residues, filter_list=filter_list,
                                          kmer_output=kmer_output)
            if first:
                labels += string_vectorize(return_labels=True, kmer=kmer, start=start, end=end, map_function=map_function, residues=residues, filter_list=filter_list)
                if features_output_format == "simple":
                    output_features(format="matrix", output_filename=features_output_filebase, labels=labels)
                    
            first = False
            i += 1

            #print(features)

            if features_output_format == "simple":
                # we'll output as we go- this is especially good for very large input files
                output_features(feature_sets=[features,], format="matrix", output_filename=features_output_filebase, write_mode="a")
            
            # we'll output sieve patterns as we go to provide a record
            if features_output_format in ("sieve", "both"):
                output_features(feature_sets=[features,], format="sieve", output_filename=features_output_filebase,
                                example_index=example_index, write_mode="a")

            if features_output_format != "simple":
                # we only keep this if we are not already dumping each line to an output file
                feature_sets.append(features)

    handle.close()

    if features_output_format != "simple":
        output_features(feature_sets=feature_sets, format=features_output_format, output_filename=features_output_filebase, example_index=example_index, labels=labels)
        
if __name__ == "__main__":
    optdict, infiles = process_options(OPTION_LIST)
    if optdict["parameterfile"]:
        optdict = read_parameterfile(**optdict)
    main(*infiles, **optdict)
