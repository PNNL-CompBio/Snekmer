Overview
========

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
Kmers are commonly used to represent sequences, particularly DNA or RNA sequences, as
counts vectors; in other words, kmer vectors simply store the number of times a given
kmer appears in a sequence. This representation is commonly used to rapidly compare
sequences based on the similarity of their kmer profiles and create vectorized
representation of sequences that are amenable to machine learning-based approaches.

Amino Acid Recoding (AAR)
-------------------------

Amino acid recoding (AAR) is a process by which the set of 20 canonical amino acids
is mapped to a smaller character set by grouping amino acids based on chemical properties
or similar measures. This technique effectively reduces the space of kmers (subsequences)
for sequence-to-sequence comparison. Proteins can thus be represented as simplified vectors
that encode structural information relevant for functional differentiation. The combined
AAR-kmer approach was successfully applied to develop machine learning (ML) models capable
of classifying sets of proteins that are functionally related but with little sequence 
similarity that would thus not be detected via conventional similarity assessment techniques [1].

.. _alphabets:

Alphabets 
`````````

Snekmer comes with 6 different recoding schemes, or "alphabets", by which AAR can occur.
The available options, including the option to skip recoding, are listed below.

=============  ===============  ======  ===========================================================================================
 Alphabet No.   Alphabet Name    Size                                         Description  
=============  ===============  ======  ===========================================================================================
 ``0``         ``hydro``          2      2-value hydrophobicity alphabet :footcite:p:`Arnold2009`
-------------  ---------------  ------  -------------------------------------------------------------------------------------------
 ``1``         ``standard``       7      "Standard" reduction alphabet :footcite:p:`Arnold2009`
-------------  ---------------  ------  -------------------------------------------------------------------------------------------
 ``2``         ``solvacc``        3      Solvent accessibility alphabet :footcite:p:`Arnold2009`
-------------  ---------------  ------  -------------------------------------------------------------------------------------------
 ``3``         ``hydrocharge``    3      2-value hydrophobicity with charged residues as a third category; by @biodataganache
-------------  ---------------  ------  -------------------------------------------------------------------------------------------
 ``4``         ``hydrostruct``    3      2-value hydrophobicity with structural-breakers as a third category; by @biodataganache
-------------  ---------------  ------  -------------------------------------------------------------------------------------------
 ``5``         ``miqs``           10     MIQS alphabet :footcite:p:`Yamada2014`
-------------  ---------------  ------  -------------------------------------------------------------------------------------------
 *n/a*         ``None``           20     No reduced alphabet
=============  ===============  ======  ===========================================================================================

Citations
:::::::::

.. footbibliography::
