#!/usr/bin/env python3
# Snekmer Model/Cluster/Search demo (Python, cross-platform).
# Run from: repo_root/resources/model_cluster_search_demo
# Steps:
#   1) Reset demo workspace
#   2) Copy query FASTAs and config
#   3) Run snekmer cluster
#   4) Run snekmer model
#   5) Collect model artifacts (optional convenience)
#   6) Run snekmer search

import sys
import shutil
import subprocess
from pathlib import Path

HERE = Path(__file__).resolve().parent
WORK = HERE
INPUT = WORK / "input"
OUTPUT = WORK / "output"

DEMO_INPUTS = HERE.parent / "demo_sequences" / "model_cluster_search_inputs"
CONFIG_SRC = HERE.parent / "config.yaml"

def show_dir(path: Path):
    print(f"Directory: {path}")
    if path.is_dir():
        for p in sorted(path.iterdir()):
            print(p.name)
    else:
        print("(missing)")
    print()

def run(cmd):
    print("+ " + " ".join(cmd))
    subprocess.run(cmd, check=True, cwd=str(WORK))

def cp(src: Path, dst: Path):
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(str(src), str(dst))

def main():
    print("Step 1: Reset demo workspace")
    for d in (OUTPUT, INPUT, WORK / ".snakemake"):
        shutil.rmtree(d, ignore_errors=True)
    INPUT.mkdir(parents=True, exist_ok=True)

    print("Step 2: Copy query FASTAs and config")
    for name in ("nxrA.faa", "nirS.faa", "TIGR03149.faa"):
        src = DEMO_INPUTS / name
        if not src.exists():
            print(f"ERROR: Missing input: {src}")
            sys.exit(2)
        cp(src, INPUT / name)

    if not CONFIG_SRC.exists():
        print(f"ERROR: Missing config: {CONFIG_SRC}")
        sys.exit(2)
    cp(CONFIG_SRC, WORK / "config.yaml")

    # show sources and workspace state
    show_dir(DEMO_INPUTS)
    show_dir(INPUT)
    show_dir(WORK)

    print("Step 3: Run snekmer cluster")
    run(["snekmer", "cluster", "--configfile=./config.yaml"])
    show_dir(OUTPUT)
    show_dir(OUTPUT / "cluster")  # if produced by your config
    show_dir(OUTPUT / "kmerize")  # if produced by your config

    print("Step 4: Run snekmer model")
    run(["snekmer", "model", "--configfile=./config.yaml"])
    show_dir(OUTPUT / "model")
    show_dir(OUTPUT / "scoring")

    print("Step 5: Collect model artifacts (optional)")
    example = OUTPUT / "example-model"
    example.mkdir(parents=True, exist_ok=True)
    for subdir, globpat in (
        ("model", "*model"),
        ("kmerize", "*kmers"),
        ("scoring", "*scorer"),
    ):
        srcdir = OUTPUT / subdir
        if srcdir.is_dir():
            for f in srcdir.glob(globpat):
                if f.is_file():
                    shutil.copy2(str(f), str(example / f.name))
    show_dir(example)

    print("Step 6: Run snekmer search")
    run(["snekmer", "search", "--configfile=./config.yaml"])
    show_dir(OUTPUT / "search")
    show_dir(OUTPUT)

    print("Done. Artifacts in: output/ (and output/example-model/).")

if __name__ == "__main__":
    main()
