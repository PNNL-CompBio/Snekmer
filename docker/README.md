README.md
---------

# Snekmer Docker Container Usage  

This directory contains instructions for a Docker setup that builds Snekmer using Python 3.10 and a virtual environment (venv) as well an example run command for using Snekmer.

## Build  

Default branch (main):
```
    docker build -t snekmer .
```
Choose a different branch at build time:
```
    docker build --build-arg SNEKMER_BRANCH=pre_paper_updates -t snekmer .
```
## Run  

Use your local directory as the working directory inside the container:

```
    docker run --rm -it -v "$PWD":/work snekmer --help
```
  
Run a command (example: learn):
```
    docker run --rm -it -v "$PWD":/work snekmer learn
```
In this case, you must have your directory setup like you normally would for running snekmer learn. For example, you will need `inputs` and `annotations` directories and a `config` file.  
   
## Notes 
- The container’s default working directory is /work.  
- The -it flags keep colored Snakemake logs.  