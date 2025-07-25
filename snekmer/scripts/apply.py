# ---------------------------------------------------------
# Imports
# ---------------------------------------------------------

import snekmer as skm
from snekmer.apply import KmerCompare

# ---------------------------------------------------------
# Files and Parameters
# ---------------------------------------------------------
config = snakemake.config

# ---------------------------------------------------------
# Run script
# ---------------------------------------------------------

weight_top = config["learn_apply"].get("weight_top", 0.5)
weight_distance = config["learn_apply"].get("weight_distance", 0.5)
save_apply_associations = config["learn_apply"]["save_apply_associations"]
        
apply = KmerCompare(
    snakemake.input.compare_associations,
    snakemake.input.data,
    snakemake.input.confidence_associations,
    snakemake.input.decoy_stats,
    snakemake.output.seq_ann,
    snakemake.output.kmer_summary,
    snakemake.params.selection_type,
    snakemake.params.threshold_type,
)
apply.execute_all(weight_top, weight_distance, save_apply_associations)
