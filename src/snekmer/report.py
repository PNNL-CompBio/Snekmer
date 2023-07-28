"""report: Snekmer report generation code.

author: @abbyjerger, @christinehc

"""
from jinja2 import Environment, PackageLoader
import os

env = Environment(loader=PackageLoader("snekmer"), auto_reload=False)

TEMPLATES = {
    "cluster": "cluster_template.html",
    "model": "model_template.html",
    "search": "search_template.html",
}


def correct_rel_path(filepath: str, out_dir: str = "output") -> str:
    """Correct relative file paths for Jinja.

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
    # force paths to use OS separator
    filepath = filepath.replace("/", os.sep)

    filepath = filepath.split(os.sep)
    out_dir = out_dir.split(os.sep)
    for d in out_dir:
        try:
            filepath.remove(d)
        except ValueError:  # if output dirs not found, skip
            pass
    return os.sep.join(filepath)


# cluster can use just this function since the image outputs are fixed at 3
# model and search also eventually call this function
def create_report(template_vars, template: str, report_file_name: str):
    """Create Snekmer reports.

    Parameters
    ----------
    template_vars : dict
        Variables defined in Jinja template
    template : str
        Name of template ("cluster", "model", or "search")
    report_file_name : str
        /path/to/report_file.html

    Returns
    -------
    None
        Creates file and exits.

    """
    # look in this file folder for the templates
    template = env.get_template(TEMPLATES[template])
    html = template.render(template_vars)
    with open(report_file_name, "w") as f:
        f.write(html)


# model
def create_report_many_images(
    image_path: str,
    table_path: str,
    rep_vars: dict,
    template: str,
    report_file_name: str,
):
    """Create report for Snekmer model mode.

    Parameters
    ----------
    image_path : str
        /path/to/images/
    table_path : str
        /path/to/tables/
    rep_vars : dict
        Variables defined in Jinja template
    template : str
        Name of template ("cluster", "model", or "search")
    report_file_name : str
        /path/to/report_file.html

    Returns
    -------
    None
        Creates file and exits.

    """
    image_list, table_list = list(), list()
    # collate images
    for root, dirs, files in os.walk(image_path):
        for file in files:
            if file.endswith(".png"):
                image_list.append(os.path.join(root, file))

    # collate tables
    for root, dirs, files in os.walk(table_path):
        for file in files:
            if file.endswith(".csv"):
                table_list.append(os.path.join(root, file))
    rep_vars["images"] = [correct_rel_path(f) for f in image_list]
    rep_vars["tables"] = [correct_rel_path(f) for f in table_list]

    create_report(rep_vars, template, report_file_name)


def create_report_many_csvs(
    path: str, rep_vars: dict, template: str, report_file_name: str
):
    # create lists with directory names and file names
    file_list = []
    dirs_list = []
    for root, dirs, files in os.walk(path):
        for d in dirs:
            dirs_list.append(d)
        for file in files:
            if file.endswith(".csv"):
                file_list.append(os.path.join(root, file))

    # combine both lists into one, in the correct order to print in the jinja template
    ordered = []
    for i in dirs_list:
        ordered.append(i)
        for j in file_list:
            if os.path.split(os.path.split(j)[0])[1] == i:
                ordered.append(j)

    rep_vars['dirs_list'] = dirs_list
    rep_vars['ordered'] = [correct_rel_path(f) for f in ordered]
    create_report(rep_vars, template, report_file_name)
