"""unzip.smk: Module for zipped file handling.

author: @christinehc
"""

# imports
import snekmer as skm

# define rules
rule unzip:
    input:
        lambda wildcards: join("input", f"{wildcards.uz}.{uz_map[wildcards.uz]}.gz")
        # lambda wildcards: unzip(join("input", ))
    output:
        join("input", "{uz}.{uzext}")
    params:
        outdir=join("input", 'zipped')
    run:

# shell:
#     "mkdir {params.outdir} && cp {input} {params.outdir} && gunzip -c {input} > {output}"
        if wildcards.sample.endswith('.fastq'):
            shell("echo gzip {input}")
            shell("echo mv {input}.gz {params.outdir}")
        else:
            shell("mv {input} {params.outdir}")
