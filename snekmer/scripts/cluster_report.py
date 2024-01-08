# imports
from os.path import dirname, join
import snekmer as skm

# check for figures
if str(snakemake.config["cluster"]["cluster_plots"]) == "True":
    fig_params = {
        "image1_name": "PCA Explained Variance",
        "image1_path": skm.report.correct_rel_path(
            join(snakemake.input.figs, "pca_explained_variance_curve.png")
        ),
        "image2_name": "Clusters (UMAP)",
        "image2_path": skm.report.correct_rel_path(
            join(snakemake.input.figs, "umap.png")
        ),
        "image3_name": "Clusters (t-SNE)",
        "image3_path": skm.report.correct_rel_path(
            join(snakemake.input.figs, "tsne.png")
        ),
    }
else:
    fig_params = {
        "image1_name": "",
        "image1_path": None,
        "image2_name": "",
        "image2_path": None,
        "image3_name": "",
        "image3_path": None,
    }

# cluster
cluster_vars = {
    "page_title": "Snekmer Cluster Report",
    "title": "Snekmer Cluster Results",
    "text": (
"Snekmer clustering results are linked below. "
        "If `cluster_plots` are enabled in the config, "
        "they will be shown below."
    ),
    "dir": dirname(skm.report.correct_rel_path(snakemake.input.table)),
    "table": skm.report.correct_rel_path(snakemake.input.table),
    **fig_params,
}

skm.report.create_report(cluster_vars, "cluster", snakemake.output[0])