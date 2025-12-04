#!/usr/bin/env bash
# Snekmer Learn→Apply demo.
# Run from repo_root/learn_apply_demo.
# This script provides directory snapshots along the way to help the user follow along.
# Steps:
#   1) Clear previous runs
#   2) Copy inputs, annotations, and config
#   3) Run snekmer learn
#   4) Prepare apply handoff files
#   5) Run snekmer apply

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

show_key_subdirs() {
  show_dir "learn"
  show_dir "learn/annotations"
  show_dir "learn/input"
  show_dir "learn/output"
  show_dir "learn/output/learn"
  show_dir "learn/output/eval_conf"

  show_dir "apply"
  show_dir "apply/annotations"
  show_dir "apply/input"
  show_dir "apply/counts"
  show_dir "apply/confidence"
  show_dir "apply/stats"
  show_dir "apply/output"
}

echo "Step 1: Clear previous runs"
rm -rf learn apply .snakemake


echo "Step 2: Create directories, copy inputs, annotations, config.yaml to learn and apply"
mkdir -p learn/input apply/input
mkdir -p learn/annotations apply/annotations

cp ../config.yaml learn/config.yaml
cp ../config.yaml apply/config.yaml

cp ../demo_sequences/learn_apply_inputs/annotations/TIGRFAMs_annotation.ann learn/annotations/
cp ../demo_sequences/learn_apply_inputs/annotations/TIGRFAMs_annotation.ann apply/annotations/

cp ../demo_sequences/learn_apply_inputs/learn/training_sequences_*.fasta learn/input/
cp ../demo_sequences/learn_apply_inputs/apply/test_sequences_1.fasta apply/input/

# concise snapshot
show_dir "learn"
show_dir "learn/annotations"
show_dir "learn/input"
show_dir "apply"
show_dir "apply/annotations"
show_dir "apply/input"

echo "Step 3: Run learn"
echo "Running snekmer learn using 10 fasta files as training input..."
cd learn
snekmer learn --configfile=./config.yaml
cd ..

# show just the new learn outputs
show_dir "learn/output"
show_dir "learn/output/learn"
show_dir "learn/output/eval_conf"

echo "Step 4: Prepare apply handoff files"
mkdir -p apply/counts apply/confidence apply/stats
cp learn/output/learn/kmer_counts_total.csv            apply/counts/
cp learn/output/eval_conf/global_confidence_scores.csv apply/confidence/
cp learn/output/eval_conf/family_summary_stats.csv     apply/stats/

# concise snapshot of apply inputs + handoff
show_dir "apply/counts"
show_dir "apply/confidence"
show_dir "apply/stats"

echo "Step 5: Run apply"
echo "Running snekmer apply using 1 fasta files as test input..."
cd apply
snekmer apply --configfile=./config.yaml
cd ..

# final concise snapshot of apply outputs
show_dir "apply/output"

echo "Done. Learn outputs: learn/output ; Apply outputs: apply/output"

echo "\n\nSample of final apply results:"

head apply/snekmer_results.csv
