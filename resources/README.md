# Resources

## YAMLs

The example YAML files included are:

- `config.yaml`: Configuration file for running Snekmer
- `clust.yaml`: (optional) Cluster configuration file for
  deploying Snekmer on a high-performance computing (HPC) cluster
  - Note: The parameter `account` must be filled in with the
    relevant account details.

## Tutorials

The `tutorial` folder contains a Jupyter notebook outlining certain aspects
of Snekmer, as well as a folder containing example files with code to run
a demo example.

### Demonstration Example

To run the demo example, run the
[demo example code](https://github.com/PNNL-CompBio/Snekmer/tree/main/resources/tutorial/demo_example):

```bash
bash run_demo.sh
```