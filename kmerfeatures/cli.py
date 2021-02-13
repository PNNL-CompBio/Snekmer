"""cli: Command line interface for Kmer pipeline."""

# imports
import argparse

from multiprocessing import cpu_count
from pkg_resources import resource_filename
from snakemake import snakemake, parse_config

# define options
MAP_FN_DESC = ["Hydrophobic/hydrophilic", "Standard 7",
               "Solvent accessibility", "LoveHateCharge",
               "LoveHateBadstruct"]
MAP_FNS = {f'reduced_alphabet_{n}': MAP_FN_DESC[n] for n in range(5)}
FEAT_OUT_FMTS = {'simple': 'Simple tab-delimited format',
                 'gist': 'Input for gist',
                 'sieve': 'Input for sieve'}

def main():
    parser = argparse.ArgumentParser(description='A program to generate features and run SIEVE models on input sequences')

    # kmer options -- remove in favor of config.yaml?
    # parser.add_argument('-v', '--verbose', action='store_false', help='verbose output')
    # parser.add_argument('-p', '--parameterfile', metavar='PATH', help='filename of a parameter file to use')
    # parser.add_argument('-f', '--fastafile', metavar='PATH', help='FASTA-format file containing protein sequences')
    # parser.add_argument('-k', '--kmer', metavar='N', type=int, default=3, help='kmer length for features')
    # parser.add_argument('-s', '--start', metavar='N', type=int, help='first residue of sequences to start generating features')
    # parser.add_argument('-n', '--end', metavar='N', type=int, help='last residue of sequences to start generating features')
    # parser.add_argument('-d', '--nucleotide', action='store_true', help='process nucleotide sequences')
    # parser.add_argument('-M', '--map_function', metavar='FUNCTION', choices=MAP_FNS, help=('mapping function for reduced amino acid alphabets:' + str(MAP_FNS)))
    # parser.add_argument('-r', '--randomize_alphabet', action='store_true', help='randomize alphabet used - but using the same number of categories and distribution as specificed by the map_function')
    # parser.add_argument('-R', '--min_rep_thresh', type=float, default=1, help='minimum number of sequences to include feature for prefiltering. 0>R<1 is percentage of input sequences')
    # parser.add_argument('-e', '--example_indexfile', metavar='FPATH', help='file containing identifiers of positive examples for sieve output records or gist .class files')
    # parser.add_argument('-m', '--features_output_format', metavar='OPTION', type=str, default='simple', choices=FEAT_OUT_FMTS.keys(), help=('format for outputting feature sets' + str(FEAT_OUT_FMTS)))
    # parser.add_argument('-o', '--features_output_filebase', metavar='BASENAME', help='filename base (no suffix) for output features')
    # parser.add_argument('-D', '--filter_duplicates', action='store_false', help='filter out duplicate identifiers from the fasta files and/or blast output')
    # parser.add_argument('-Q', '--output_shuffled_sequences', metavar='N', type=int, help='output shuffled sequences in a fasta file for each input sequence')
    # parser.add_argument('-F', '--feature_set', help='specify which features to include in the model')
    # parser.add_argument('-K', '--kmer_output', action='store_false', help='verbose output of kmer positions for each sequence [under development]')
    # parser.add_argument('-w', '--walk', action='store_true', help='perform a kmer walk on the input fasta file to get an idea of the kmer representation')

    # snakemake options
    parser.add_argument('--dryrun', action='store_true', help='perform a dry run')
    parser.add_argument('--configfile', metavar='PATH', default='config.yaml', help='path to yaml configuration file')
    parser.add_argument('--config', nargs="*", metavar="KEY=VALUE",
                        help=("Set or overwrite values in the"
                              " workflow config object. The workflow"
                              " config object is accessible as"
                              " variable config inside the workflow."
                              " Default values can be set by"
                              " providing a JSON file."
                              )
                        )
    parser.add_argument('--unlock', action='store_true', help='unlock directory')
    parser.add_argument('--touch', action='store_true', help='touch output files only')
    parser.add_argument('--latency', metavar='N', type=int, default=3, help='specify filesystem latency (seconds)')
    parser.add_argument('--cores', metavar='N', type=int, default=cpu_count(), help='number of cores used for execution (local execution only)')
    parser.add_argument('--count', metavar='N', type=int, help='number of files to process (limits DAG size)')
    parser.add_argument('--countstart', metavar='IDX', type=int, default=0, help='starting file index (for use with --count)')

    # cluster options
    clust = parser.add_argument_group('cluster arguments')
    clust.add_argument('--cluster', metavar='PATH', help='path to cluster execution yaml configuration file')
    clust.add_argument('--jobs', metavar='N', type=int, default=1000, help='number of simultaneous jobs to submit to a slurm queue')

    # parse args
    args = parser.parse_args()
    config = parse_config(args)

    # start/stop config
    if args.count is not None:
        config = {'start': args.countstart,
                  'stop': args.countstart + args.count,
                  **config}
    else:
        config = {**config}

    # cluster config
    if args.cluster is not None:
        cluster = "sbatch -A {cluster.account} -N {cluster.nodes} -t {cluster.time} -J {cluster.name} --ntasks-per-node {cluster.ntasks} -p {cluster.partition}"
    else:
        cluster = None

    snakemake(resource_filename('kmerfeatures', 'Snakefile'),
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
              touch=args.touch,
              latency_wait=args.latency)


if __name__ == '__main__':
    main()
