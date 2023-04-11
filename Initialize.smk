import os
import sys
from warnings import warn
from py.helpers import *
configfile: "config.yaml"

rule_all = [
    os.path.join(config["OUT_DIR"], "sample_info.tsv"),  # data 
    expand(
        os.path.join(
            config["SNP_DIR"],
            "1kGP_high_coverage_Illumina.{chr}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
        ),
        chr=get_chromosomes(config["CHROMOSOMES"])
    ),
    expand(
        os.path.join(
            config["SNP_DIR"],
            "1kGP_high_coverage_Illumina.{chr}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz.tbi"
        ),
        chr=get_chromosomes(config["CHROMOSOMES"]),
    ),
    expand(
        os.path.join(
            config["SNP_TR_PANEL_DIR"],
            "{chr}_TR_SNP_merged.vcf.gz"
        ),
        chr=get_chromosomes(config["CHROMOSOMES"])
    ),
    expand(
        os.path.join(
            config["SNP_TR_PANEL_DIR"],
            "{chr}_TR_SNP_merged.vcf.gz.tbi"
        ),
        chr=get_chromosomes(config["CHROMOSOMES"])
    ),
    os.path.join(config["OUT_DIR"], "sample_info.tsv"),
    os.path.join(config["OUT_DIR"], "sample_info_common.tsv")
]

rule all: input: rule_all

rule download_snps:
    output: 
        vcf = os.path.join(config["SNP_DIR"], "1kGP_high_coverage_Illumina.{chr}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"),
        tbi = os.path.join(config["SNP_DIR"], "1kGP_high_coverage_Illumina.{chr}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz.tbi")
    params:
        chr = "{chr}",
    threads: 1
    resources: mem_mb=1000
    log: "logs/download.snp.{chr}.log"
    shell:
        """
        mkdir -p {params.out_dir}/1kG_SNP
        wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.{params.chr}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz \
            -O {output.vcf}
        wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.{params.chr}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz.tbi \
            -O {output.tbi}
        """


rule download_SNP_TR_panel:
    output:
        vcf = os.path.join(config["SNP_TR_PANEL_DIR"], "{chr}_TR_SNP_merged.vcf.gz"),
        tbi = os.path.join(config["SNP_TR_PANEL_DIR"], "{chr}_TR_SNP_merged.vcf.gz.tbi")
    params:
        chr="{chr}"
    threads: 1
    resources: mem_mb=1000
    log: "logs/download.tr.{chr}.log"
    shell:
        """
        wget https://ensemble-tr.s3.us-east-2.amazonaws.com/phased-split/{params.chr}_final_SNP_merged.vcf.gz \
            -O {output.vcf}
        wget https://ensemble-tr.s3.us-east-2.amazonaws.com/phased-split/{params.chr}_final_SNP_merged.vcf.gz.csi \
            -O {output.tbi}
        """


rule make_sample_info_file:
    input: os.path.join(config["SNP_DIR"], "1kGP_high_coverage_Illumina.chr21.filtered.SNV_INDEL_SV_phased_panel.vcf.gz")
    output: 
        samp_info = os.path.join(config["OUT_DIR"], "sample_info.tsv"),
        common_info = os.path.join(config["OUT_DIR"], "sample_info_common.tsv")
    params:
        data_dir = config["DATA_DIR"],
        ignore = config["IGNORE_SAMPLES"]
    threads: 1
    resources:
        mem_mb = 50
    log: "logs/mk.info.log"
    shell:
        """
        ls -1 {params.data_dir}/*.bam \
            | tr '\\n' '\\0' \
            | xargs -0 -n 1 basename \
            | tr '.' '\\t' \
            > {output.samp_info}
        
        >{output.common_info}
        while read SID SEQTYPE HAPLOTYPE SORTED BAM; do
            if [[ $(bcftools query -l {input}) == *"${{SID}}"* ]] && [[ "{params.ignore}" != *"${{SID}}"* ]]; then
                echo "${{SID}}\t${{SEQTYPE}}\t${{HAPLOTYPE}}\t${{SORTED}}\t${{BAM}}" >> {output.common_info}
            fi
        done <{output.samp_info}
        """
