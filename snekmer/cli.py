"""cli: Command line interface for Snekmer.

author: @christinehc, @biodataganache

"""
# imports
import argparse

from multiprocessing import cpu_count
from os.path import join
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
    parser = {}

    # main parser
    parser["main"] = argparse.ArgumentParser(
        description="Snekmer: A tool for kmer-based sequence analysis using amino acid reduction (AAR)"
    )
    parser["main"].add_argument(
        "-v",
        "--version",
        action="version",
        version=__version__,
        help="print version and exit",
    )

    # snakemake options
    parser["smk"] = argparse.ArgumentParser(add_help=False)
    parser["smk"].add_argument(
        "--dryrun", "-n", action="store_true", help="perform a dry run"
    )
    parser["smk"].add_argument(
        "--configfile",
        metavar="PATH",
        default="config.yaml",
        help="path to yaml configuration file",
    )
    parser["smk"].add_argument(
        "-C",
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
    parser["smk"].add_argument("--unlock", action="store_true", help="unlock directory")
    parser["smk"].add_argument(
        "-U",
        "--until",
        metavar="TARGET",
        nargs="+",
        help="run pipeline until reaching the target rule or files",
    )
    parser["smk"].add_argument(
        "--keepgoing",
        "--keep-going",
        "-k",
        action="store_true",
        default=False,
        help="go on with independent jobs if a job fails",
    )
    parser["smk"].add_argument(
        "--latency",
        "-w",
        "--output-wait",
        "--latency-wait",
        metavar="SECONDS",
        type=int,
        default=30,
        help="wait time, in seconds, for output file creation (default 30)",
    )
    parser["smk"].add_argument(
        "--touch", action="store_true", help="touch output files only"
    )
    parser["smk"].add_argument(
        "-c",
        "--cores",
        metavar="N",
        type=int,
        default=cpu_count(),
        help="number of cores used for execution (local execution only)",
    )
    parser["smk"].add_argument(
        "--count",
        metavar="N",
        type=int,
        help="number of files to process (limits DAG size)",
    )
    parser["smk"].add_argument(
        "--countstart",
        metavar="IDX",
        type=int,
        default=0,
        help="starting file index (for use with --count)",
    )
    parser["smk"].add_argument(
        "--verbose",
        action="store_true",
        default=False,
        help="show additional debug output (default False)",
    )

    # clust execution options
    parser["clust"] = parser["smk"].add_argument_group("cluster execution arguments")
    parser["clust"].add_argument(
        "--clust",
        metavar="PATH",
        help="path to cluster execution yaml configuration file",
    )
    parser["clust"].add_argument(
        "-j",
        "--jobs",
        metavar="N",
        type=int,
        default=1000,
        help="number of simultaneous jobs to submit to a slurm queue",
    )

    # create subparsers for each operation mode
    # parser.add_argument("mode", choices=["cluster", "model", "search"])
    parser["subparsers"] = parser["main"].add_subparsers(
        title="mode", description="Snekmer mode", dest="mode"
    )

    # subparsers
    parser["cluster"] = parser["subparsers"].add_parser(
        "cluster",
        description="Apply unsupervised clustering via Snekmer",
        parents=[parser["smk"]],
    )
    parser["model"] = parser["subparsers"].add_parser(
        "model",
        description="Train supervised models via Snekmer",
        parents=[parser["smk"]],
    )
    parser["search"] = parser["subparsers"].add_parser(
        "search",
        description="Search sequences against pre-existing models via Snekmer",
        parents=[parser["smk"]],
    )

    # parse args
    args = parser["main"].parse_args()
    config = parse_config(args)

    # start/stop config
    if args.count is not None:
        config = {
            "start": args.countstart,
            "stop": args.countstart + args.count,
            **config,
        }
    else:
        config = config

    # cluster config
    if args.clust is not None:
        cluster = "sbatch -A {cluster.account} -N {cluster.nodes} -t {cluster.time} -J {cluster.name} --ntasks-per-node {cluster.ntasks} -p {cluster.partition}"
    else:
        cluster = None

    # parse operation mode
    if args.mode == "cluster":
        snakemake(
            resource_filename("snekmer", join("rules", "cluster.smk")),
            configfiles=[args.configfile],
            config=config,
            cluster_config=args.clust,
            cluster=cluster,
            keepgoing=args.keepgoing,
            force_incomplete=True,
            cores=args.cores,
            nodes=args.jobs,
            dryrun=args.dryrun,
            unlock=args.unlock,
            until=args.until,
            touch=args.touch,
            latency_wait=args.latency,
            verbose=args.verbose,
        )

    elif args.mode == "model":
        snakemake(
            resource_filename("snekmer", join("rules", "model.smk")),
            configfiles=[args.configfile],
            config=config,
            cluster_config=args.clust,
            cluster=cluster,
            keepgoing=args.keepgoing,
            force_incomplete=True,
            cores=args.cores,
            nodes=args.jobs,
            dryrun=args.dryrun,
            unlock=args.unlock,
            until=args.until,
            touch=args.touch,
            latency_wait=args.latency,
            verbose=args.verbose,
        )

    elif args.mode == "search":
        snakemake(
            resource_filename("snekmer", join("rules", "search.smk")),
            configfiles=[args.configfile],
            config=config,
            cluster_config=args.clust,
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
            verbose=args.verbose,
        )

    else:
        parser["main"].print_help()


if __name__ == "__main__":
    main()
