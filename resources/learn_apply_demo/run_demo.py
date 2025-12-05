#!/usr/bin/env python3
# Snekmer Learn→Apply demo (Python). Run from repo_root/resources/learn_apply_demo.

import os, sys, shutil, subprocess, glob
from pathlib import Path

HERE = Path(__file__).resolve().parent
LEARN = HERE / "learn"
APPLY = HERE / "apply"

CONFIG_SRC = HERE.parent / "config.yaml"
ANN_SRC = HERE.parent / "demo_sequences" / "learn_apply_inputs" / "annotations" / "TIGRFAMs_annotation.ann"
LEARN_SRC_DIR = HERE.parent / "demo_sequences" / "learn_apply_inputs" / "learn"
APPLY_SRC_FILE = HERE.parent / "demo_sequences" / "learn_apply_inputs" / "apply" / "test_sequences_1.fasta"

def show_dir(p: Path):
    print(f"Directory: {p}")
    if p.is_dir():
        for n in sorted(os.listdir(p)):
            print(n)
    else:
        print("(missing)")
    print()

# Step 1
print("Step 1: Clear previous runs")
for d in (LEARN, APPLY, HERE / ".snakemake"):
    shutil.rmtree(d, ignore_errors=True)

# Step 2
print("\nStep 2: Create directories, copy inputs, annotations, config.yaml to learn and apply")
(LEARN / "input").mkdir(parents=True, exist_ok=True)
(LEARN / "annotations").mkdir(parents=True, exist_ok=True)
(APPLY / "input").mkdir(parents=True, exist_ok=True)
(APPLY / "annotations").mkdir(parents=True, exist_ok=True)

# configs
shutil.copy2(CONFIG_SRC, LEARN / "config.yaml")
shutil.copy2(CONFIG_SRC, APPLY / "config.yaml")

# annotations (ensure in learn/annotations)
shutil.copy2(ANN_SRC, LEARN / "annotations" / ANN_SRC.name)
shutil.copy2(ANN_SRC, APPLY / "annotations" / ANN_SRC.name)

# learn inputs
learn_glob = sorted(glob.glob(str(LEARN_SRC_DIR / "training_sequences_*.fasta")))
if not learn_glob:
    raise FileNotFoundError(f"No files matched: {LEARN_SRC_DIR}/training_sequences_*.fasta")
for src in learn_glob:
    shutil.copy2(src, LEARN / "input" / Path(src).name)

# apply input
shutil.copy2(APPLY_SRC_FILE, APPLY / "input" / APPLY_SRC_FILE.name)

# concise snapshot
show_dir(LEARN)
show_dir(LEARN / "annotations")
show_dir(LEARN / "input")
show_dir(APPLY)
show_dir(APPLY / "annotations")
show_dir(APPLY / "input")

# Step 3
print("Step 3: Run learn")
print("Running snekmer learn using 10 fasta files as training input...")
subprocess.run(["snekmer", "learn", "--configfile=./config.yaml"], check=True, cwd=str(LEARN))

show_dir(LEARN / "output")
show_dir(LEARN / "output" / "learn")
show_dir(LEARN / "output" / "eval_conf")

# Step 4
print("Step 4: Prepare apply handoff files")
(APPLY / "counts").mkdir(parents=True, exist_ok=True)
(APPLY / "confidence").mkdir(parents=True, exist_ok=True)
(APPLY / "stats").mkdir(parents=True, exist_ok=True)

shutil.copy2(LEARN / "output" / "learn" / "kmer_counts_total.csv",            APPLY / "counts" / "kmer_counts_total.csv")
shutil.copy2(LEARN / "output" / "eval_conf" / "global_confidence_scores.csv", APPLY / "confidence" / "global_confidence_scores.csv")
shutil.copy2(LEARN / "output" / "eval_conf" / "family_summary_stats.csv",     APPLY / "stats" / "family_summary_stats.csv")

show_dir(APPLY / "counts")
show_dir(APPLY / "confidence")
show_dir(APPLY / "stats")

# Step 5
print("Step 5: Run apply")
print("Running snekmer apply using 1 fasta files as test input...")
subprocess.run(["snekmer", "apply", "--configfile=./config.yaml"], check=True, cwd=str(APPLY))

show_dir(APPLY / "output")

print("Done. Learn outputs: learn/output ; Apply outputs: apply/output")
print("\n\nSample of final apply results:")
try:
    with (APPLY / "snekmer_results.csv").open("r", encoding="utf-8", errors="replace") as fh:
        for i, line in enumerate(fh):
            if i >= 10: break
            sys.stdout.write(line)
except FileNotFoundError:
    print("(missing) apply/snekmer_results.csv")
