configfile: "config.yaml"
from py.helpers import *

sample_dict = get_sample_dict(config, init=False)

rule all:
    input: expand(os.path.join(config["DATA_DIR"], "pileup", "{sample}.pileup"), sample=sample_dict.keys())

rule get_pileup:
    input: os.path.join(config["DATA_DIR"], "short", "{sample}.short{ext}")
    output: os.path.join(config["DATA_DIR"], "pileup", "{sample}.pileup")
    resources:
        mem_mb=4000
    threads: 1
    params:
        ref=config["REF_FASTA"]
    conda: "envs/sam.yaml"
    shell: "samtools mpileup -f {params.ref} {input} > {output}"
