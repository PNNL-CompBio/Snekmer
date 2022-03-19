"""cli: Command line interface for Snekmer.

author: @christinehc, @biodataganache

"""
# imports
import argparse

from multiprocessing import cpu_count
from pkg_resources import resource_filename
from snakemake import snakemake, parse_config
from snekmer import __version__

# define options
MAP_FN_DESC = [
    "Hydrophobic/hydrophilic",
    "Standard 7",
    "Solvent accessibility",
    "LoveHateCharge",
    "LoveHateBadstruct",
]
MAP_FNS = {f"reduced_alphabet_{n}": MAP_FN_DESC[n] for n in range(5)}
FEAT_OUT_FMTS = {
    "simple": "Simple tab-delimited format",
    "gist": "Input for gist",
    "sieve": "Input for sieve",
}


def main():
    parser = argparse.ArgumentParser(
        description="Snekmer: A program to generate features and run SIEVE models on input sequences"
    )
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version=__version__,
        help="print version and exit",
    )
    # parser.add_argument('--mode', choices=['model', 'cluster'], help='select operation mode')

    # snakemake options
    parser.add_argument("--dryrun", action="store_true", help="perform a dry run")
    parser.add_argument(
        "--configfile",
        metavar="PATH",
        default="config.yaml",
        help="path to yaml configuration file",
    )
    parser.add_argument(
        "--config",
        nargs="*",
        metavar="KEY=VALUE",
        help=(
            "Set or overwrite values in the"
            " workflow config object. The workflow"
            " config object is accessible as"
            " variable config inside the workflow."
            " Default values can be set by"
            " providing a JSON file."
        ),
    )
    parser.add_argument("--unlock", action="store_true", help="unlock directory")
    parser.add_argument(
        "--until",
        metavar="TARGET",
        nargs="+",
        help="run pipeline until reaching the target rule or files",
    )
    parser.add_argument(
        "--latency",
        "--latency-wait",
        metavar="SECONDS",
        type=int,
        default=30,
        help="wait time, in seconds, for output file creation (default 30)",
    )
    parser.add_argument("--touch", action="store_true", help="touch output files only")
    parser.add_argument(
        "--cores",
        metavar="N",
        type=int,
        default=cpu_count(),
        help="number of cores used for execution (local execution only)",
    )
    parser.add_argument(
        "--count",
        metavar="N",
        type=int,
        help="number of files to process (limits DAG size)",
    )
    parser.add_argument(
        "--countstart",
        metavar="IDX",
        type=int,
        default=0,
        help="starting file index (for use with --count)",
    )

    # cluster options
    clust = parser.add_argument_group("cluster arguments")
    clust.add_argument(
        "--cluster",
        metavar="PATH",
        help="path to cluster execution yaml configuration file",
    )
    clust.add_argument(
        "--jobs",
        metavar="N",
        type=int,
        default=1000,
        help="number of simultaneous jobs to submit to a slurm queue",
    )

    # create subparsers for both operation modes
    parser.add_argument("mode", choices=["model", "cluster", "search"])

    # parse args
    args = parser.parse_args()
    config = parse_config(args)

    # start/stop config
    if args.count is not None:
        config = {
            "start": args.countstart,
            "stop": args.countstart + args.count,
            "mode": args.mode,
            **config,
        }
    else:
        config = {"mode": args.mode, **config}

    # cluster config
    if args.cluster is not None:
        cluster = "sbatch -A {cluster.account} -N {cluster.nodes} -t {cluster.time} -J {cluster.name} --ntasks-per-node {cluster.ntasks} -p {cluster.partition}"
    else:
        cluster = None

    # parse operation mode
    if args.mode == "model":
        snakemake(
            resource_filename("snekmer", "rules/model.smk"),
            configfiles=[args.configfile],
            config=config,
            cluster_config=args.cluster,
            cluster=cluster,
            keepgoing=True,
            force_incomplete=True,
            cores=args.cores,
            nodes=args.jobs,
            dryrun=args.dryrun,
            unlock=args.unlock,
            until=args.until,
            touch=args.touch,
            latency_wait=args.latency,
        )

    elif args.mode == "cluster":
        snakemake(
            resource_filename("snekmer", "rules/cluster.smk"),
            configfiles=[args.configfile],
            config=config,
            cluster_config=args.cluster,
            cluster=cluster,
            keepgoing=True,
            force_incomplete=True,
            cores=args.cores,
            nodes=args.jobs,
            dryrun=args.dryrun,
            unlock=args.unlock,
            until=args.until,
            touch=args.touch,
            latency_wait=args.latency,
        )

    elif args.mode == "search":
        snakemake(
            resource_filename("snekmer", "rules/search.smk"),
            configfiles=[args.configfile],
            config=config,
            cluster_config=args.cluster,
            cluster=cluster,
            keepgoing=True,
            force_incomplete=True,
            cores=args.cores,
            nodes=args.jobs,
            dryrun=args.dryrun,
            unlock=args.unlock,
            until=args.until,
            touch=args.touch,
            latency_wait=args.latency,
        )


if __name__ == "__main__":
    main()
