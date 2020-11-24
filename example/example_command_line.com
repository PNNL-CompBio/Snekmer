python3 KmerFeatures.py -f example/example.fasta -k 14 -M reduced_alphabet_0 -m simple -o example/example_features_1 -F example/model_features.txt
