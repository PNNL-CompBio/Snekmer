# ---------------------------------------------------------
# Imports
# ---------------------------------------------------------

import snekmer as skm
from snekmer.apply import KmerCompare

# ---------------------------------------------------------
# Helper for config
# ---------------------------------------------------------

def _extract_learn_apply_cfg(config):
    """
    Safely extract the learn_apply sub-config.

    Supports two cases:
      - config is a full Snakemake-style config dict with a "learn_apply" key
      - config is already just the learn_apply dict
    """
    if config is None:
        return {}
    # Snakemake-style: top-level config with "learn_apply"
    if isinstance(config, dict) and "learn_apply" in config:
        la_cfg = config.get("learn_apply") or {}
    else:
        # Assume caller passed only the learn_apply section
        la_cfg = config
    return la_cfg or {}


def run_apply(
    compare_associations,
    data,
    confidence_associations,
    decoy_stats,
    seq_ann_out,
    kmer_summary_out,
    selection_type,
    threshold_type,
    config=None,
):
    """
    Run the KmerCompare apply step.

    Parameters
    ----------
    compare_associations : str
        Path to the kmer-association totals CSV.
    data : str
        Path to the input NPZ (sequences/k-mers to annotate).
    confidence_associations : str
        Path to the global confidence CSV (Difference → confidence).
    decoy_stats : str
        Path to reverse/decoy family stats CSV.
    seq_ann_out : str
        Output path for per-sequence annotations.
    kmer_summary_out : str
        Output path for the summary k-mer counts.
    selection_type : str
        Selection method name (e.g., "best", "hybrid", etc.).
    threshold_type : str
        Threshold column name (e.g., "Median", "Percentile90", or "None").
    config : dict or None
        Full Snakemake config or just the "learn_apply" section.
    """
    la_cfg = _extract_learn_apply_cfg(config)

    weight_top = la_cfg.get("weight_top", 0.5)
    weight_distance = la_cfg.get("weight_distance", 0.5)
    # Original behavior: required key; here we default to False only if missing
    save_apply_associations = la_cfg.get("save_apply_associations", False)

    apply = KmerCompare(
        compare_associations,
        data,
        confidence_associations,
        decoy_stats,
        seq_ann_out,
        kmer_summary_out,
        selection_type,
        threshold_type,
    )
    apply.execute_all(weight_top, weight_distance, save_apply_associations)


# ---------------------------------------------------------
# Snakemake entry point (kept for pipeline compatibility)
# ---------------------------------------------------------

if "snakemake" in globals():
    cfg = snakemake.config
    run_apply(
        compare_associations=snakemake.input.compare_associations,
        data=snakemake.input.data,
        confidence_associations=snakemake.input.confidence_associations,
        decoy_stats=snakemake.input.decoy_stats,
        seq_ann_out=snakemake.output.seq_ann,
        kmer_summary_out=snakemake.output.kmer_summary,
        selection_type=snakemake.params.selection_type,
        threshold_type=snakemake.params.threshold_type,
        config=cfg,
    )
