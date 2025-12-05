#!/usr/bin/env python3
# snekmer/scripts/apply_concat.py
# Concatenate per-sample k-mer summaries into a single CSV.
# ---------------------------------------------------------
# Imports
# ---------------------------------------------------------
from typing import Sequence
import os

def concat_csv(inputs: Sequence[str], output: str) -> None:
    """
    Concatenate CSVs by writing the header from the first file,
    then appending all subsequent rows from every file (skipping each header).

    Parameters
    ----------
    inputs : Sequence[str]
        Ordered list of CSV paths to concatenate.
    output : str
        Output CSV path to write.
    """
    if not inputs:
        os.makedirs(os.path.dirname(output) or ".", exist_ok=True)
        open(output, "w", encoding="utf-8").close()
        return

    os.makedirs(os.path.dirname(output) or ".", exist_ok=True)

    first = inputs[0]
    with open(output, "w", encoding="utf-8", newline="") as w:
        # Write header from first file
        with open(first, "r", encoding="utf-8", newline="") as r0:
            header = r0.readline()
            if header:
                w.write(header)

        # Append data rows from every file, skipping the header line each time
        for f in inputs:
            with open(f, "r", encoding="utf-8", newline="") as r:
                _ = r.readline()  # skip header
                for line in r:
                    w.write(line)

# ---------------------------------------------------------
# Snakemake entry point (pipeline-only)
# ---------------------------------------------------------

if "snakemake" in globals():
    inp = [str(p) for p in snakemake.input]
    out = str(snakemake.output[0])
    concat_csv(inp, out)
