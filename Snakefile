import os
import sys
from warnings import warn
from py.helpers import *
configfile: "config.yaml"

rule_all = [
    expand(
        os.path.join(config["OUT_DIR"],
            "tricolor", 
            "sensor", 
            "{sid}_{seqtype}_sensor"
        ), zip,
        sid=get_samp_id(os.path.join(config["OUT_DIR"], "sample_info_common.tsv")),
        seqtype=get_seqtype(os.path.join(config["OUT_DIR"], "sample_info_common.tsv"))
    ),
    expand(
        os.path.join(config["OUT_DIR"],
            "tricolor", 
            "refer", 
            "{sid}_{seqtype}_refer"
        ), zip,
        sid=get_samp_id(os.path.join(config["OUT_DIR"], "sample_info_common.tsv")),
        seqtype=get_seqtype(os.path.join(config["OUT_DIR"], "sample_info_common.tsv"))
    ),
    expand(
        os.path.join(config["OUT_DIR"],
            "straglr", 
            "{sid}_{seqtype}_straglr",
            "{sid}_{seqtype}_hp1_straglr.tsv"
        ), zip,
        sid=get_samp_id(os.path.join(config["OUT_DIR"], "sample_info_common.tsv")),
        seqtype=get_seqtype(os.path.join(config["OUT_DIR"], "sample_info_common.tsv"))
    ),
    expand(
        os.path.join(config["OUT_DIR"],
            "straglr", 
            "{sid}_{seqtype}_straglr",
            "{sid}_{seqtype}_hp2_straglr.tsv"
        ), zip,
        sid=get_samp_id(os.path.join(config["OUT_DIR"], "sample_info_common.tsv")),
        seqtype=get_seqtype(os.path.join(config["OUT_DIR"], "sample_info_common.tsv"))
    ),
    expand(
        os.path.join(config["OUT_DIR"],
            "imputed_TR_SNPs",
            "{chr}"
        ),
        chr=get_chromosomes(config["CHROMOSOMES"]),
    )

]

print(get_seqtype(os.path.join(config["OUT_DIR"], "sample_info_common.tsv")))

rule all: input: rule_all


rule tricolor:
    output: 
        sensor = os.path.join(config["OUT_DIR"], "tricolor", "sensor", "{sid}_{seqtype}_sensor"),
        refer = os.path.join(config["OUT_DIR"], "tricolor", "refer", "{sid}_{seqtype}_refer")
    input:
        hp1 = os.path.join(config["DATA_DIR"], "{sid}.{seqtype}.h1.sorted.bam"),
        hp2 = os.path.join(config["DATA_DIR"], "{sid}.{seqtype}.h2.sorted.bam")
    params: 
        data_dir = os.path.join(config["OUT_DIR"]),
        ref_fasta = config["REF_FASTA"],
        sid = "{sid}",
        seqtype = "{seqtype}"
    threads: 1
    resources: mem_mb=2000
    log: "logs/tricolor_{sid}_{seqtype}.log"
    conda: "stromboli.yml"
    shell:
        """
        mkdir -p {params.data_dir}/tricolor/sensor/{params.sid}_{params.seqtype}_sensor
        mkdir -p {params.data_dir}/tricolor/refer
        mkdir -p {params.data_dir}/tricolor/sage

        TRiCoLOR SENSoR -bam {input.hp1} {input.hp2} -o {output.sensor}
        TRiCoLOR REFER -g {params.ref_fasta} -bam {input.hp1} {input.hp2} -bed {output.sensor}/*.bed -o {output.refer}
        """


rule straglr:
    input:
        hp1 = os.path.join(config["DATA_DIR"], "{sid}.{seqtype}.h1.sorted.bam"),
        hp2 = os.path.join(config["DATA_DIR"], "{sid}.{seqtype}.h2.sorted.bam")
    output: 
        hp1 = os.path.join(config["OUT_DIR"], "straglr", "{sid}_{seqtype}_straglr", "{sid}_{seqtype}_hp1_straglr.tsv"),
        hp2 = os.path.join(config["OUT_DIR"], "straglr", "{sid}_{seqtype}_straglr", "{sid}_{seqtype}_hp2_straglr.tsv")
    params: 
        out_dir = os.path.join(config["OUT_DIR"], "straglr", "{sid}_{seqtype}_straglr"),
        ref_fasta = config["REF_FASTA"]
    threads: 1
    resources: mem_mb=2000
    log: "logs/straglr.{sid}.{seqtype}.log"
    conda: "stromboli.yml"
    shell:
        """
        mkdir -p {params.out_dir}/strglr_temp
        python straglr.py \
            {input.hp1} \
            {params.ref_fasta} \
            "$(basename {input.hp1} .bam)" \
            --tmpdir {params.out_dir}/strglr_temp
        
        python straglr.py \
            {input.hp2} \
            {params.ref_fasta} \
            "$(basename {input.hp2} .bam)" \
            --tmpdir {params.out_dir}/strglr_temp
        """


rule impute:
    input: 
        snp = os.path.join(config["SNP_DIR"], "1kGP_high_coverage_Illumina.{chr}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"),
        tr_snp = os.path.join(config["SNP_TR_PANEL_DIR"], "{chr}_TR_SNP_merged.vcf.gz"),
        common_samples = os.path.join(config["OUT_DIR"], "sample_info_common.tsv")
    output: os.path.join(config["OUT_DIR"], "imputed_TR_SNPs", "{chr}")
    params:
        out_dir = os.path.join(config["OUT_DIR"], "straglr", "{chr}_straglr"),
        beagle_jar = config["BEAGLE_JAR_PATH"],
        chr = "{chr}"

    threads: 4 
    resources: mem_mb=8000
    log: "logs/impute.{chr}.log"
    conda: "stromboli.yml"
    shell:
        """
        mkdir -p {params.out_dir}/imputed_TRs
        java -Xmx8g -jar {params.beagle_jar} \
            gt=$(bcftools view -S {input.common_samples} {input.snp})  \
            ref=$(bcftools view -S {input.common_samples} {input.tr_snp}) \
            out={output}
        """