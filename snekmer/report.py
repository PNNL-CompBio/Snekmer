from jinja2 import Environment, PackageLoader
import os

env = Environment(loader=PackageLoader("snekmer"), auto_reload=False)

TEMPLATES = {
    "cluster": "cluster_template.html",
    "model": "model_template.html",
    "search": "search_template.html",
}


def correct_rel_path(filepath: str, out_dir: str = "output") -> str:
    """Correct relative file paths for Jinja

    Parameters
    ----------
    filepath : str
        /path/to/file
    out_dir : str, optional
        /path/to/output, by default "output"

    Returns
    -------
    str
        /path/to/file
    """
    filepath = filepath.split(os.sep)
    out_dir = out_dir.split(os.sep)
    for d in out_dir:
        filepath.remove(d)
    return os.sep.join(filepath)


# cluster can use just this function since the image outputs are fixed at 3
# model and search also eventually call this function
def create_report(template_vars, template, report_file_name):
    """Create Snekmer reports.

    Parameters
    ----------
    template_vars : dict
        _description_
    template : str
        Name of template ("cluster", "model", or "search")
    report_file_name : _type_
        _description_

    """
    # look in this file folder for the templates
    template = env.get_template(TEMPLATES[template])
    html = template.render(template_vars)
    with open(report_file_name, "w") as f:
        f.write(html)


# model or search
def create_report_many_images(path, rep_vars, template, report_file_name):
    """Create report for Snekmer model or search modes.

    Parameters
    ----------
    path : _type_
        _description_
    rep_vars : _type_
        _description_
    template : _type_
        _description_
    report_file_name : _type_
        _description_
    """
    filelist = []
    for root, dirs, files in os.walk(path):
        for file in files:
            if file.endswith(".png"):
                filelist.append(os.path.join(root, file))
    rep_vars["images"] = filelist
    create_report(rep_vars, template, report_file_name)


# cluster
cluster_vars = dict(
    page_title="Snekmer Cluster Report",
    title="Snekmer Cluster Results",
    image1_name="PCA",
    image1_path="/path/to/pca_explained_variance.png",
    image2_name="t-SNE",
    image2_path="/path/to/tsne.png",
    image3_name="UMAP",
    image3_path="/path/to/umap.png",
    text="These are the three figure outputs produced by the Snekmer Cluster command.",
)
cluster_template = "cluster_template.html"
cluster_report = "Snekmer_Cluster_Report.html"

# create_report(cluster_vars, cluster_template, cluster_report)

# model
model_vars = dict(
    page_title="Snekmer Model Report",
    title="Snekmer Model Results",
    text="These are the figure outputs produced by the Snekmer Model command.",
)
model_template = "model_template.html"
model_report = "Snekmer_Model_Report.html"

# create_report_many_images("output/model/figures", model_vars,
#                           model_template, model_report)

# search
search_vars = dict(
    page_title="Snekmer Search Report",
    title="Snekmer Search Results",
    text="These are the outputs file locations from Snekmer Search.",
    dir="/path/to/search/",
)
search_template = "search_template.html"
search_report = "Snekmer_Search_Report.html"

# create_report(search_vars, search_template, search_report)
