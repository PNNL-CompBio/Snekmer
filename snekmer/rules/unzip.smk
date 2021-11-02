"""unzip.smk: Module for zipped file handling.

author: @christinehc
"""

# imports
from os.path import join
import snekmer as skm

# define rules
rule unzip:
    input:
        join("input", "{unzipped}.gz")
        # lambda wildcards: join("input", f"{wildcards.uz}.{uz_map[wildcards.uz]}.gz")
        # lambda wildcards: unzip(join("input", ))
    output:
        join("input", "{unzipped}")
    params:
        outdir=join("input", 'zipped')
    shell:
        "mkdir {params.outdir} && cp {input} {params.outdir} && gunzip -c {input} > {output}"
    # run:
        # if wildcards.sample.endswith('.fastq'):
        #     shell("echo gzip {input}")
        #     shell("echo mv {input}.gz {params.outdir}")
        # else:
        #     shell("mv {input} {params.outdir}")
