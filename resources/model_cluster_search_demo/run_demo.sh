#!/usr/bin/env bash
# Snekmer Model/Cluster/Search demo.
# Run from repo_root/model_cluster_search_demo.
# Steps:
#   1) Reset demo workspace
#   2) Copy query FASTAs and config
#   3) Run snekmer cluster
#   4) Run snekmer model
#   5) Collect model artifacts (optional convenience)
#   6) Run snekmer search

set -e

show_dir() {
  local d="$1"
  echo "Directory: $d"
  if [[ -d "$d" ]]; then
    ls -1 "$d" || true
  else
    echo "(missing)"
  fi
  echo
}

echo "Step 1: Reset demo workspace"
rm -rf output input .snakemake
mkdir -p input

echo "Step 2: Copy query FASTAs and config"
cp ../demo_sequences/model_cluster_search_inputs/nxrA.faa input/
cp ../demo_sequences/model_cluster_search_inputs/nirS.faa input/
cp ../demo_sequences/model_cluster_search_inputs/TIGR03149.faa input/

cp ../config.yaml ./config.yaml

ls -1 ../demo_sequences/model_cluster_search_inputs/ || true
echo
show_dir "input"
show_dir "."

echo "Step 3: Run snekmer cluster"
snekmer cluster --configfile=./config.yaml
show_dir "output"
show_dir "output/cluster"      # if produced by your config
show_dir "output/kmerize"      # if produced by your config

echo "Step 4: Run snekmer model"
snekmer model --configfile=./config.yaml
show_dir "output/model"
show_dir "output/scoring"

echo "Step 5: Collect model artifacts (optional)"
mkdir -p output/example-model
# copy (not move) so subsequent steps still find expected paths
cp -f output/model/*model   output/example-model/ 2>/dev/null || true
cp -f output/kmerize/*kmers output/example-model/ 2>/dev/null || true
cp -f output/scoring/*scorer output/example-model/ 2>/dev/null || true
show_dir "output/example-model"

echo "Step 6: Run snekmer search"
snekmer search --configfile=./config.yaml
show_dir "output/search"     
show_dir "output"

echo "Done. Artifacts in: output/ (and output/example-model/)."