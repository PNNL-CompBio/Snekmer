# !/bin/bash

# reset files
mv input/zipped/*.gz input/
rm -rf output

# run snekmer model on examples using provided config.yaml
snekmer model --configfile=../../motif_config.yaml

# move the model files to bespoke directories
mkdir output/example-model
mv output/model/*model output/example-model/
mv output/kmerize/*kmers output/example-model/
mv output/scoring/*scorer output/example-model/

# run snekmer motif on examples using provided config.yaml
snekmer motif --configfile=../../motif_config.yaml
