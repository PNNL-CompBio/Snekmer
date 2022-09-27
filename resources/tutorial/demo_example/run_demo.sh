#!/bin/bash
# run snekmer cluster on examples using provided config.yaml
snekmer cluster --configfile=../../config.yaml

# run snekmer model on examples using provided config.yaml
snekmer model --configfile=../../config.yaml

# move the model files to bespoke directories
mkdir output/example-model
mv output/model/*model output/example-model/
mv output/kmerize/*kmers output/example-model/
mv output/scoring/*scorer output/example-model/

# run snekmer search on examples using provided config.yaml
snekmer search --configfile=../../config.yaml

