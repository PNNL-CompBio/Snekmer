# KmerPipeline
Pipeline to apply encoded Kmer analysis to protein sequences

Input: fasta protein sequences
Output: clusters of similar proteins

Evaluation output: assessment of how well the clusters of similar proteins represent functions

Steps envisioned:
1. Kmer feature generation (using KmerFeatures.py)
2. Similarity analysis - based on kmer vectors - Jaccard similarity?
3. Clustering of similarity graph to communities - MCL? Other?

Evaluation pipeline:
1. Start with well-annotated protein sequences
2. Parse annotations from protein sequences into usable form (if neessary)
3. Run pipeline above with some parameters
4. Assess how sensitive the method is - that is how likely is it that a members of a cluster have the annotation that is the most prevalent for that cluster
5. Assess the specificity - that is how likely is it that an annotation maps to a single cluster