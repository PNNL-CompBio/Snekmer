"""Shared pytest fixtures for Snekmer tests."""
import gzip
import textwrap
from pathlib import Path

import pytest


# ---------------------------------------------------------------------------
# Minimal FASTA content helpers
# ---------------------------------------------------------------------------

_FASTA_UNIPROT = textwrap.dedent("""\
    >sp|FAM001|SEQ1_HUMAN Protein one OS=Homo sapiens
    MKTIIALSYIFCLVFA
    >sp|FAM001|SEQ2_HUMAN Protein two OS=Homo sapiens
    ACDEFGHIKLMNPQRSTVWY
    >sp|FAM002|SEQ3_MOUSE Protein three OS=Mus musculus
    MKTIIALSYIFCLVFAACDE
""")

_FASTA_NO_PIPES = textwrap.dedent("""\
    >SEQ1 plain header no pipes
    MKTIIALSYIFCLVFA
    >SEQ2 another plain header
    ACDEFGHIKLMNPQRSTVWY
""")

_ANN_CONTENT = textwrap.dedent("""\
    id\tfamily
    FAM001\tFAM001
    FAM002\tFAM002
""")


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def fasta_file(tmp_path) -> Path:
    """Single FASTA file with UniProt-style headers (two families, three seqs)."""
    p = tmp_path / "seqs.fasta"
    p.write_text(_FASTA_UNIPROT)
    return p


@pytest.fixture
def fasta_dir(tmp_path) -> Path:
    """Directory of three FASTA files, each with one sequence."""
    d = tmp_path / "fastas"
    d.mkdir()
    sequences = [
        ("file1.fasta", ">sp|FAM001|A Human\nMKTIIALSYIFCLVFA\n"),
        ("file2.faa",   ">sp|FAM001|B Human\nACDEFGHIKLMNPQRS\n"),
        ("file3.fna",   ">sp|FAM002|C Mouse\nMKTIIALSYIFCLVFAACDE\n"),
    ]
    for name, content in sequences:
        (d / name).write_text(content)
    return d


@pytest.fixture
def fasta_gz_file(tmp_path) -> Path:
    """Gzip-compressed FASTA file."""
    p = tmp_path / "seqs.fasta.gz"
    with gzip.open(p, "wt", encoding="utf-8") as fh:
        fh.write(_FASTA_UNIPROT)
    return p


@pytest.fixture
def fasta_no_pipes(tmp_path) -> Path:
    """FASTA file whose headers have no pipe characters."""
    p = tmp_path / "nopipes.fasta"
    p.write_text(_FASTA_NO_PIPES)
    return p


@pytest.fixture
def ann_file(tmp_path) -> Path:
    """Minimal .ann file (tab-separated id/family)."""
    p = tmp_path / "annotations.ann"
    p.write_text(_ANN_CONTENT)
    return p


@pytest.fixture
def fasta_empty(tmp_path) -> Path:
    """Empty file (no sequences)."""
    p = tmp_path / "empty.fasta"
    p.write_text("")
    return p


@pytest.fixture
def fasta_dna(tmp_path) -> Path:
    """FASTA file whose sequences are DNA (ATCGN >90%)."""
    p = tmp_path / "dna.fasta"
    p.write_text(
        ">seq1 a nucleotide sequence\nATCGATCGATCGATCGATCG\n"
        ">seq2 another nucleotide sequence\nGCTAGCTAGCTAGCTAGCTA\n"
    )
    return p


@pytest.fixture
def ann_wrong_cols(tmp_path) -> Path:
    """Annotation file with wrong column names (no 'family')."""
    p = tmp_path / "bad_cols.ann"
    p.write_text("id\tannotation\nFAM001\tsomefam\n")
    return p


@pytest.fixture
def ann_comma(tmp_path) -> Path:
    """Comma-separated annotation file (common mistake)."""
    p = tmp_path / "comma.ann"
    p.write_text("id,family\nFAM001,FAM001\nFAM002,FAM002\n")
    return p


@pytest.fixture
def ann_one_family(tmp_path) -> Path:
    """Annotation file with only one distinct family value."""
    p = tmp_path / "one_family.ann"
    p.write_text("id\tfamily\nFAM001\tonlyfam\n")
    return p


@pytest.fixture
def ann_no_data(tmp_path) -> Path:
    """Annotation file with header only, no data rows."""
    p = tmp_path / "no_data.ann"
    p.write_text("id\tfamily\n")
    return p


@pytest.fixture
def learn_output_dir(tmp_path) -> Path:
    """Simulated learn output directory with the three handoff files."""
    learn_dir = tmp_path / "learn"
    counts_dir = learn_dir / "apply_inputs" / "counts"
    stats_dir = learn_dir / "apply_inputs" / "stats"
    conf_dir = learn_dir / "apply_inputs" / "confidence"
    for d in (counts_dir, stats_dir, conf_dir):
        d.mkdir(parents=True)
    (counts_dir / "kmer_counts_total.csv").write_text("family,k1\nFAM001,1\n")
    (stats_dir / "family_summary_stats.csv").write_text("family,Median\nFAM001,0.5\n")
    (conf_dir / "global_confidence_scores.csv").write_text("score\n0.9\n")
    return learn_dir
