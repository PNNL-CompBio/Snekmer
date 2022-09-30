from jinja2 import Environment, FileSystemLoader
import os


# cluster can use just this function since the image outputs are fixed at 3
# model and search also eventually call this function
def create_report(template_vars, template_name, report_file_name):
    # look in this file folder for the templates
    env = Environment(loader=FileSystemLoader('templates'))
    template = env.get_template(template_name)
    html = template.render(template_vars)
    with open(report_file_name, 'w') as f:
        f.write(html)


# model or search
def create_report_many_images(path, rep_vars, template_name, report_file_name):
    filelist = []
    for root, dirs, files in os.walk(path):
        for file in files:
            if file.endswith(".png"):
                filelist.append(os.path.join(root, file))
    rep_vars["images"] = filelist
    create_report(rep_vars, template_name, report_file_name)


# cluster
cluster_vars = dict(page_title='Snekmer Cluster Report',
                    title='Snekmer Cluster Results',
                    image1_name='PCA',
                    image1_path='C:\\Users\\jerg881\\SnekmerTestingCluster\\output\\cluster\\figures'
                                '\\pca_explained_variance_curve.png',
                    image2_name='t-SNE',
                    image2_path='C:\\Users\\jerg881\\SnekmerTestingCluster\\output\\cluster\\figures\\tsne.png',
                    image3_name='UMAP',
                    image3_path='C:\\Users\\jerg881\\SnekmerTestingCluster\\output\\cluster\\figures\\umap.png',
                    text='These are the three figure outputs produced by the Snekmer Cluster command.')
cluster_template = 'cluster_template.html'
cluster_report = 'Snekmer_Cluster_Report.html'

create_report(cluster_vars, cluster_template, cluster_report)

# model
model_vars = dict(page_title='Snekmer Model Report',
                  title='Snekmer Model Results',
                  text='These are the figure outputs produced by the Snekmer Model command.')
model_template = 'model_template.html'
model_report = 'Snekmer_Model_Report.html'

create_report_many_images("C:\\Users\\jerg881\\SnekmerTestingModel\\output\\model\\figures", model_vars,
                          model_template, model_report)

# search
search_vars = dict(page_title='Snekmer Search Report',
                   title='Snekmer Search Results',
                   text='These are the outputs file locations from Snekmer Search.',
                   dir='C:\\Users\\jerg881\\SnekmerTestingSearch\\output\\search')
search_template = 'search_template.html'
search_report = 'Snekmer_Search_Report.html'

create_report(search_vars, search_template, search_report)
