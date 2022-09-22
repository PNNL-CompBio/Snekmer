#!/bin/bash

# run snekmer model on examples using provided config.yaml
snekmer model --configfile=../../config.yaml

# run snekmer search on examples using provided config.yaml
snekmer search --configfile=../../config.yaml

# run snekmer cluster on examples using provided config.yaml
snekmer cluster --configfile=../../config.yaml
