This directory has an example protein fasta file and the output from running
the original KmerFeatures.py code (from the repo biodataganache/SIEVE-Ub) on it.
I've included the command that was used to generate the file:
python3 KmerFeatures.py -f example/example.fasta -k 14 -M reduced_alphabet_0 -m simple -o example/example_features_1 -F example/model_features.txt
