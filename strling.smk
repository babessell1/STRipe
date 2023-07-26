import os
import sys
from warnings import warn
from py.helpers import *
configfile: "config.yaml"

rule_all = [
    expand(
        os.path.join(config["OUT_DIR"], "strling", "{sid}.final-genotype.txt"),
        zip,
        sid=get_samp_id(os.path.join(config["OUT_DIR"], "sample_info_common.tsv")),
        seqtype=get_seqtype(os.path.join(config["OUT_DIR"], "sample_info_common.tsv"))
    ),
    expand(
        os.path.join(config["SHORT_READS_DIR"], "{sid}.final.cram.crai"),
        zip,
        sid=get_samp_id(os.path.join(config["OUT_DIR"], "sample_info_common.tsv"))
    )
]


rule all: input: rule_all

# Rule to index the input cram file
rule index_cram:
    input:
        os.path.join(config["SHORT_READS_DIR"], "{sid}.final.cram")
    output:
        os.path.join(config["SHORT_READS_DIR"], "{sid}.final.cram.crai")
    conda: "envs/sam.yaml"
    shell:
        """
        samtools index -@ 2 {input}
        """

# Rule to run strling
rule strling:
    input:
        cram=os.path.join(config["SHORT_READS_DIR"], "{sid}.final.cram"),
        index=os.path.join(config["SHORT_READS_DIR"], "{sid}.final.cram.crai")
    output:
        os.path.join(config["OUT_DIR"], "strling", "{sid}.final-genotype.txt")
    params:
        out_dir = config["OUT_DIR"],
        sid = "{sid}",
        seqtype = "{seqtype}",
        strling = config["STRLING_PATH"],
        ref_fasta = config["REF_FASTA"]
    threads: 1
    resources: mem_mb = 4000
    shell:
        """
        bname=$(basename "{input.cram}" .cram)
        {params.strling} extract -f "{params.ref_fasta}" "{input.cram}" "output/${{bname}}.bin"
        {params.strling} call --output-prefix "{params.out_dir}/strling//${{bname}}" -f "{params.ref_fasta}" "{input.cram}" "output/${{bname}}.bin"
        """