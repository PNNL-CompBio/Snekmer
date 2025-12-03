"""cli: Command line interface for Snekmer.

author: @christinehc, @biodataganache, @snakemake

"""
# imports
import argparse
import os
import sys
import json
import shlex
import subprocess

from multiprocessing import cpu_count
from importlib.resources import files, as_file
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


def get_argument_parser():
    parser = {}

    # main parser
    parser["main"] = argparse.ArgumentParser(
        description=(
            "Snekmer: A tool for kmer-based sequence analysis"
            " using amino acid reduction (AAR)."
        )
    )
    parser["main"].add_argument(
        "-v",
        "--version",
        action="version",
        version=__version__,
        help="Print version and exit.",
    )

    # Snakemake options
    # NOTE: Arguments and documentation borrowed from Snakemake:
    # https://github.com/snakemake/snakemake/blob/main/snakemake/__init__.py
    parser["smk"] = argparse.ArgumentParser(add_help=False)
    parser["smk"].add_argument(
        "--dry-run",
        "--dryrun",
        "-n",
        dest="dryrun",
        action="store_true",
        help="Do not execute anything, and display what would be done. "
        "If you have a very large workflow, use --dry-run --quiet to just "
        "print a summary of the DAG of jobs.",
    )
    parser["smk"].add_argument(
        "--configfile",
        nargs="+",
        # default="config.yaml",
        metavar="PATH",
        help=(
            "Specify or overwrite the config file of the workflow (see the docs). "
            "Values specified in JSON or YAML format are available in the global config "
            "dictionary inside the workflow. Multiple files overwrite each other in "
            "the given order. Thereby missing keys in previous config files are extended by "
            "following configfiles. Note that this order also includes a config file defined "
            "in the workflow definition itself (which will come first)."
        ),
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
        "--until",
        "-U",
        nargs="+",
        metavar="TARGET",
        help=(
            "Runs the pipeline until it reaches the specified rules or "
            "files. Only runs jobs that are dependencies of the specified "
            "rule or files, does not run sibling DAGs. "
        ),
    )
    parser["smk"].add_argument(
        "--keepgoing",
        "--keep-going",
        "-k",
        action="store_true",
        help="Go on with independent jobs if a job fails.",
    )
    parser["smk"].add_argument(
        "--latency",
        "-w",
        "--output-wait",
        "--latency-wait",
        type=int,
        default=30,
        metavar="SECONDS",
        help="Wait given seconds if an output file of a job is not present after "
        "the job finished. This helps if your filesystem "
        "suffers from latency (default 30).",
    )
    parser["smk"].add_argument(
        "--touch",
        "-t",
        action="store_true",
        help=(
            "Touch output files (mark them up to date without really "
            "changing them) instead of running their commands. This is "
            "used to pretend that the rules were executed, in order to "
            "fool future invocations of snakemake. Fails if a file does "
            "not yet exist. Note that this will only touch files that would "
            "otherwise be recreated by Snakemake (e.g. because their input "
            "files are newer). For enforcing a touch, combine this with "
            "--force, --forceall, or --forcerun. Note however that you loose "
            "the provenance information when the files have been created in "
            "realitiy. Hence, this should be used only as a last resort."
        ),
    )
    parser["smk"].add_argument(
        "--cores",
        "-c",
        default=cpu_count(),
        type=int,
        metavar="N",
        help=(
            "Use at most N CPU cores/jobs in parallel. "
            "If N is omitted or 'all', the limit is set to the number of "
            "available CPU cores. "
            "In case of cluster/cloud execution, this argument sets the maximum number "
            "of cores requested from the cluster or cloud scheduler. (See "
            "https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#"
            "resources-remote-execution for more info)"
            "This number is available to rules via workflow.cores."
        ),
    )
    parser["smk"].add_argument(
        "--count",
        metavar="N",
        type=int,
        help="Number of files to process (limits DAG size).",
    )
    parser["smk"].add_argument(
        "--countstart",
        metavar="IDX",
        type=int,
        default=0,
        help="Starting file index (for use with --count).",
    )
    parser["smk"].add_argument(
        "--verbose",
        action="store_true",
        default=False,
        help="Show additional debug output (default False)",
    )
    parser["smk"].add_argument(
        "--quiet",
        "-q",
        nargs="*",
        choices=["progress", "rules", "all"],
        default=None,
        help=(
            "Do not output certain information. If used without "
            "arguments, do not output any progress or rule "
            "information. Defining 'all' results in no "
            "information being printed at all."
        ),
    )
    parser["smk"].add_argument(
        "--directory",
        "-d",
        metavar="DIR",
        action="store",
        help=(
            "Specify working directory (relative paths in "
            "the snakefile will use this as their origin)."
        ),
    )
    parser["smk"].add_argument(
        "--forcerun",
        "-R",
        nargs="*",
        metavar="TARGET",
        help=(
            "Force the re-execution or creation of the given rules or files."
            " Use this option if you changed a rule and want to have all its "
            "output in your workflow updated."
        ),
    )
    parser["smk"].add_argument(
        "--list-code-changes",
        "--lc",
        action="store_true",
        help="List all output files for which the rule body (run or shell) have "
        "changed in the Snakefile.",
    )
    parser["smk"].add_argument(
        "--list-params-changes",
        "--lp",
        action="store_true",
        help="List all output files for which the defined params have changed "
        "in the Snakefile.",
    )

    # clust execution options
    parser["clust"] = parser["smk"].add_argument_group("Cluster Execution Arguments")
    parser["clust"].add_argument(
        "--clust",
        nargs="+",
        metavar="PATH",
        help="Path to cluster execution yaml configuration file.",
    )
    parser["clust"].add_argument(
        "-j",
        "--jobs",
        metavar="N",
        type=int,
        default=1000,
        help="Number of simultaneous jobs to submit to a slurm queue.",
    )

    # create subparsers for each operation mode
    # parser.add_argument("mode", choices=["cluster", "model", "search", "learn", "apply"])
    parser["subparsers"] = parser["main"].add_subparsers(
        title="mode",
        description="Snekmer mode (cluster, model, search, learn, or apply).",
        dest="mode",
    )

    # subparsers
    parser["cluster"] = parser["subparsers"].add_parser(
        "cluster",
        description="Apply unsupervised clustering via Snekmer.",
        parents=[parser["smk"]],
    )
    parser["model"] = parser["subparsers"].add_parser(
        "model",
        description="Train supervised models via Snekmer.",
        parents=[parser["smk"]],
    )
    parser["search"] = parser["subparsers"].add_parser(
        "search",
        description="Search sequences against pre-existing models via Snekmer.",
        parents=[parser["smk"]],
    )
    parser["motif"] = parser["subparsers"].add_parser(
        "motif",
        description="Find structurally and functionally relevant motifs via Snekmer.",
        parents=[parser["smk"]]
    )
    parser["learn"] = parser["subparsers"].add_parser(
        "learn",
        description="Learn kmer-annotation associations via Snekmer",
        parents=[parser["smk"]],
    )
    parser["apply"] = parser["subparsers"].add_parser(
        "apply",
        description="Apply kmer-annotation associations via Snekmer",
        parents=[parser["smk"]],
    )
    return parser


def get_main_args():
    parser = get_argument_parser()
    return parser["main"]


# ------------------------- Snakemake v9 bridge helpers -------------------------
def _coerce_scalar(val: str):
    """Best-effort cast of KEY=VALUE overrides (bool/int/float/json/str)."""
    low = val.lower()
    if low in {"true", "false"}:
        return low == "true"
    try:
        return int(val)
    except ValueError:
        pass
    try:
        return float(val)
    except ValueError:
        pass
    try:
        return json.loads(val)
    except Exception:
        return val


def _parse_config_overrides(args):
    """Re-implement minimal parse_config behavior (Snakemake ≥8 removed API)."""
    cfg = {}
    if args.config:
        for item in args.config:
            if "=" not in item:
                raise ValueError(f"--config expects KEY=VALUE, got: {item!r}")
            k, v = item.split("=", 1)
            cfg[k] = _coerce_scalar(v)
    # start/stop config (preserve original semantics)
    if args.count is not None:
        cfg = {"start": args.countstart, "stop": args.countstart + args.count, **cfg}
    return cfg


def _quiet_cli(args):
    """Map original quiet behavior to Snakemake CLI flags."""
    if args.quiet is None:
        return []
    if len(args.quiet) == 0:
        # default case, set quiet to progress and rule
        return ["--quiet", "progress", "--quiet", "rules"]
    out = []
    for q in args.quiet:
        out.extend(["--quiet", q])
    return out


def _resolve_rulefile(relpath):
    """Resolve packaged rules/<name>.smk using importlib.resources (replaces resource_filename)."""
    res = files("snekmer") / relpath
    return as_file(res)  # context manager yielding filesystem path


def _build_snakemake_cmd(rulefile_path, args, configfiles, config, keepgoing_override=None):
    """Translate previous snakemake(...) call into a Snakemake v9 CLI invocation."""
    cmd = ["snakemake", "-s", rulefile_path, "--cores", str(args.cores)]

    # previous parameters → CLI flags
    if args.dryrun:
        cmd.append("--dry-run")
    if args.touch:
        cmd.append("--touch")
    if args.unlock:
        cmd.append("--unlock")
    if args.list_code_changes:
        cmd.append("--list-code-changes")
    if args.list_params_changes:
        cmd.append("--list-params-changes")
    if args.forcerun:
        cmd.extend(["--forcerun", *args.forcerun])
    if args.until:
        cmd.extend(["--until", *args.until])
    if args.latency is not None:
        cmd.extend(["--latency-wait", str(args.latency)])
    if args.verbose:
        cmd.append("--verbose")

    # treat incomplete as needing re-run (force_incomplete=True)
    cmd.append("--rerun-incomplete")

    # keepgoing: allow per-mode override exactly like original code
    if keepgoing_override is True or (keepgoing_override is None and args.keepgoing):
        cmd.append("--keep-going")

    # config files
    for cf in configfiles:
        cmd.extend(["--configfiles", cf])

    # inline config pairs
    if config:
        pairs = []
        for k, v in config.items():
            if isinstance(v, (dict, list)):
                pairs.append(f"{k}={json.dumps(v)}")
            else:
                pairs.append(f"{k}={v}")
        cmd.extend(["--config", *pairs])

    # jobs
    if args.jobs is not None:
        cmd.extend(["--jobs", str(args.jobs)])

    # working directory (previously workdir=args.directory)
    if args.directory:
        cmd.extend(["--directory", os.path.abspath(args.directory)])

    # quiet flags
    cmd.extend(_quiet_cli(args))

    # Original built a fixed sbatch template when --clust was provided.
    if args.clust is not None:
        cluster = (
            "sbatch -A {cluster.account} -N {cluster.nodes} -t {cluster.time} "
            "-J {cluster.name} --ntasks-per-node {cluster.ntasks} -p {cluster.partition}"
        )
        cmd.extend(
            [
                "--executor",
                "cluster-generic",
                "--cluster-generic-submit-cmd",
                cluster,
            ]
        )
        for cc in args.clust:
            cmd.extend(["--cluster-config", cc])

    return cmd


def _run(rule_basename, args, keepgoing_override=None):
    """Helper to mirror per-mode snakemake(...) calls with minimal changes."""
    # fix configfile path
    if args.configfile is None:
        configfile = ["config.yaml"]
    else:
        configfile = args.configfile

    if (args.directory is not None) and (args.configfile is None):
        configfiles = [os.path.join(args.directory, c) for c in configfile]
    else:
        configfiles = list(map(os.path.abspath, configfile))

    # config (replacement for parse_config + start/stop handling)
    config = _parse_config_overrides(args)

    # resolve rules/<mode>.smk and invoke Snakemake CLI
    relpath = os.path.join("rules", f"{rule_basename}.smk")
    with _resolve_rulefile(relpath) as rulefile_path:
        cmd = _build_snakemake_cmd(
            str(rulefile_path), args, configfiles, config, keepgoing_override=keepgoing_override
        )
        if args.verbose:
            print(">> Running:", shlex.join(cmd), file=sys.stderr)
        proc = subprocess.run(cmd)
        return proc.returncode


def main():
    parser = get_argument_parser()

    # parse args
    args = parser["main"].parse_args()

    # set quiet output options
    if args.quiet is not None and len(args.quiet) == 0:
        # default case, set quiet to progress and rule
        args.quiet = ["progress", "rules"]

    # parse operation mode
    if args.mode == "cluster":
        rc = _run("cluster", args, keepgoing_override=args.keepgoing)

    elif args.mode == "model":
        rc = _run("model", args, keepgoing_override=args.keepgoing)

    elif args.mode == "search":

        rc = _run("search", args, keepgoing_override=True)

    elif args.mode == "motif":
        rc = _run("motif", args, keepgoing_override=args.keepgoing)

    elif args.mode == "learn":

        rc = _run("learn", args, keepgoing_override=True)

    elif args.mode == "apply":

        rc = _run("apply", args, keepgoing_override=True)

    else:
        parser["main"].print_help()
        rc = 2

    sys.exit(rc)


if __name__ == "__main__":
    main()
