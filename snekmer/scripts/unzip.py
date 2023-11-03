"""unzip.py: File unzipping handling

author: @christinehc

"""
# imports
import gzip
from os import makedirs, system
from os.path import dirname, exists


# if not exists(dirname(snakemake.output.zipped)):
#     makedirs(dirname(snakemake.output.zipped))
# system(f"mv {snakemake.input[0]} {snakemake.output.zipped}")

# unzip and save file contents
with gzip.open(snakemake.output.zipped, "rb") as openf, open(
    snakemake.output.unzipped, "wb"
) as savef:
    file_content = openf.readlines()
    for line in file_content:
        savef.write(line)

# rule perform_kmer_walk:
#     input:
#         fasta=get_fasta_files
#     output:
#         # need to fix code to properly generate an output...
#     run:
#         skm.walk.kmer_walk(input.fasta)
