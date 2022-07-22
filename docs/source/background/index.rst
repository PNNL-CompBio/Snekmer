Background
==========

Protein function annotation typically involves pairwise or multiple sequence alignment,
followed by various search tools to query for sequence similarity across pre-annotated
sequences. Snekmer aims to simplify this process by combining the kmer approach
with amino acid recoding (AAR) to map the 20 amino acids into alphabets of smaller size.
This combined approach represents protein sequences as AAR-kmer vectors, which preserve
information about the parent sequence while also reducing the sheer size of kmer vectors
based on the canonical amino acids (i.e. :math:`20^{k}`).

What is a kmer?
---------------

Kmers (or k-mers) are all of the subsequences of length *k* that comprise a sequence.

Amino Acid Recoding (AAR)
-------------------------

*to be continued...*