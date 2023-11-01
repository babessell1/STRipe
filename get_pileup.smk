import os
configfile: "config.yaml"
from py.helpers import *

sample_dict = get_sample_dict(config, init=False)

rule_all = []

rule_all.extend(
    expand(os.path.join(config["DATA_DIR"], "pileup", "{sample}.pileup"),
        zip,
        sample=get_samples(sample_dict, "short"),
        ext=get_ext(sample_dict, "short")
    )
)

rule all:
    input: rule_all


rule get_pileup:
    # input os.path.join(config["DATA_DIR"], "short", "{sample}.short{ext}") need wildcard ext
    input: os.path.join(config["DATA_DIR"], "short", "{sample}.short.cram")
    output: os.path.join(config["DATA_DIR"], "pileup", "{sample}.pileup")
    resources:
        mem_mb=4000
    threads: 1
    params:
        ref=config["REF_FASTA"]
    shell: "samtools mpileup -f {params.ref} {input.bam} > {output}"

