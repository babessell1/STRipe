import os
import sys
from warnings import warn
from py.helpers import *
configfile: "config.yaml"

rule_all = [
    os.path.join(config["OUT_DIR"], "sample_info.tsv"),  # data 
    expand(
        os.path.join(
            config["OUT_DIR"],
            "1kG_SNP",
            "1kGP_high_coverage_Illumina.{chr}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
        ),
        chr=get_chromosomes(config["CHROMOSOMES"])
    ),
    expand(
        os.path.join(
            config["OUT_DIR"],
            "1kG_SNP",
            "1kGP_high_coverage_Illumina.{chr}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz.tbi"
        ),
        chr=get_chromosomes(config["CHROMOSOMES"]),
    ),
    expand(
        os.path.join(
            config["OUT_DIR"],
            "TR_SNP_panel",
            "{chr}_TR_SNP_merged.vcf.gz"
        ),
        chr=get_chromosomes(config["CHROMOSOMES"])
    ),
    expand(
        os.path.join(
            config["OUT_DIR"],
            "TR_SNP_panel",
            "{chr}_TR_SNP_merged.vcf.gz.tbi"
        ),
        chr=get_chromosomes(config["CHROMOSOMES"])
    ),
    os.path.join(config["OUT_DIR"], "sample_info.tsv"),
    os.path.join(config["OUT_DIR"], "sample_info_common.tsv"),
    expand(
        os.path.join(config["OUT_DIR"],
            "tricolor", 
            "sensor", 
            "{sid}_{seqtype}_sensor"
        ),
        sid=get_samp_id(os.path.join(config["OUT_DIR"], "sample_info_common.tsv")),
        seqtype=get_seqtype(os.path.join(config["OUT_DIR"], "sample_info_common.tsv"))
    ),
    expand(
        os.path.join(config["OUT_DIR"],
            "tricolor", 
            "refer", 
            "{sid}_{seqtype}_refer"
        ),
        sid=get_samp_id(os.path.join(config["OUT_DIR"], "sample_info_common.tsv")),
        seqtype=get_seqtype(os.path.join(config["OUT_DIR"], "sample_info_common.tsv"))
    ),
    expand(
        os.path.join(config["OUT_DIR"],
            "straglr", 
            "{sid}_{seqtype}_straglr",
            "{sid}_{seqtype}_hp1_straglr.tsv"
        ),
        sid=get_samp_id(os.path.join(config["OUT_DIR"], "sample_info_common.tsv")),
        seqtype=get_seqtype(os.path.join(config["OUT_DIR"], "sample_info_common.tsv"))
    ),
    expand(
        os.path.join(config["OUT_DIR"],
            "straglr", 
            "{sid}_{seqtype}_straglr",
            "{sid}_{seqtype}_hp2_straglr.tsv"
        ),
        sid=get_samp_id(os.path.join(config["OUT_DIR"], "sample_info_common.tsv")),
        seqtype=get_seqtype(os.path.join(config["OUT_DIR"], "sample_info_common.tsv"))
    ),
    expand(
        os.path.join(config["OUT_DIR"],
            "imputed_TR_SNPs",
            "{sid}_{seqtype}"
        ),
        sid=get_samp_id(os.path.join(config["OUT_DIR"], "sample_info_common.tsv")),
        seqtype=get_seqtype(os.path.join(config["OUT_DIR"], "sample_info_common.tsv"))
    )

]


rule all: input: rule_all


rule download_snps:
    output: 
        vcf = os.path.join(config["OUT_DIR"], "1kG_SNP", "1kGP_high_coverage_Illumina.{chr}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"),
        tbi = os.path.join(config["OUT_DIR"], "1kG_SNP", "1kGP_high_coverage_Illumina.{chr}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz.tbi")
    params:
        chr = "{chr}",
        out_dir = os.path.join(config["OUT_DIR"], "straglr", "{sid}_{seqtype}_straglr")
    threads: 1
    resources: mem_mb=2000
    log: "./logs/dwnld.snp.log"
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
        vcf = os.path.join(config["OUT_DIR"], "TR_SNP_panel", "{chr}_TR_SNP_merged.vcf.gz"),
        tbi = os.path.join(config["OUT_DIR"], "TR_SNP_panel", "{chr}_TR_SNP_merged.vcf.gz.tbi")
    params:
        chr="{chr}"
    threads: 1
    resources: mem_mb=2000
    log: "./logs/dwnld.tr.log"
    shell:
        """
        wget https://ensemble-tr.s3.us-east-2.amazonaws.com/phased-split/{params.chr}_final_SNP_merged.vcf.gz \
            -O {output.vcf}
        wget https://ensemble-tr.s3.us-east-2.amazonaws.com/phased-split/{params.chr}_final_SNP_merged.vcf.gz.csi \
            -O {output.tbi}
        """


rule make_sample_info_file:
    input: os.path.join(config["OUT_DIR"], "snps", "1kGP_high_coverage_Illumina.chr21.filtered.SNV_INDEL_SV_phased_panel.vcf.gz")
    output: 
        samp_info = os.path.join(config["OUT_DIR"], "sample_info.tsv")
        common_info = os.path.join(config["OUT_DIR"], "sample_info_common.tsv")
    threads: 1
    resources:
        mem_mb = 50
    log: "./logs/mk.info.log"
    shell:
        """
        ls -1 {params.data_dir}/*.bam \
            | tr '\\n' '\\0' \
            | xargs -0 -n 1 basename \
            | tr '.' '\\t' \
            > {output.samp_info}
        
        >{output.common_info}
        while read SID SEQTYPE HAPLOTYPE SORTED BAM; do
            if [[ "$(vcf-query -l {input})" == *"${{SID}}"*]]; then
                echo "${SID} >> {output.common_info}
            fi
        done <{output.samp_info}
        """


rule tricolor:
    output: 
        sensor = os.path.join(config["OUT_DIR"], "tricolor", "sensor", "{sid}_{seqtype}_sensor"),
        refer = os.path.join(config["OUT_DIR"], "tricolor", "refer", "{sid}_{seqtype}_refer"),
        #app = os.path.join(config["OUT_DIR"], "tricolor", "app", "{sid}_{seqtype}_app")
    input:
        hp1 = os.path.join(config["DATA_DIR"], "{sid}.{seqtype}.h1.sorted.bam"),
        hp2 = os.path.join(config["DATA_DIR"], "{sid}.{seqtype}.h2.sorted.bam")
    params: 
        data_dir = os.path.join(config["OUT_DIR"]),
        ref_fasta = cong["REF_FASTA"],
        sid = "{sid}",
        seqtype = "{seqtype}"
    threads: 1
    resources: mem_mb=2000
    log: "./logs/tricolor.log"
    shell:
        """
        mkdir -p {params.data_dir}/tricolor/sensor/{params.sid}_{params.seqtype}_sensor
        mkdir -p {params.data_dir}/tricolor/refer
        mkdir -p {params.data_dir}/tricolor/sage
        mkdir -p {params.data_dir}/tricolor/app

        TRiCoLOR SENSoR -bam {input.hp1} {input.hp2} -o {output.sensor}
        TRiCoLOR REFER -g {params.ref_fasta} -bam {input.hp1} {input.hp2} -bed {output.sensor}/*.bed -o {output.refer}
        #TRiCoLOR ApP -g {params.ref_fasta} \
        #    -bam {output.refer}/haplotype1/TRiCoLOR.bam {output.refer}/haplotype2/TRiCoLOR.bam \
        #    -gb {output.refer}/reference/TRiCoLOR.srt.bed.gz \
        #    -h1b {output.refer}/haplotype1/TRiCoLOR.srt.bed.gz -h2b {output.refer}/haplotype2/TRiCoLOR.srt.bed.gz \
        #    -o {output.app} <REGION>
        """


rule straglr:
    input:
        hp1 = os.path.join(config["DATA_DIR"], "{sid}.{seqtype}.h1.sorted.bam"),
        hp2 = os.path.join(config["DATA_DIR"], "{sid}.{seqtype}.h2.sorted.bam")
    output: 
        hp1 = os.path.join(config["OUT_DIR"], "straglr", "{sid}_{seqtype}_straglr", "{sid}_{seqtype}_hp1_straglr.tsv")
        hp2 = os.path.join(config["OUT_DIR"], "straglr", "{sid}_{seqtype}_straglr", "{sid}_{seqtype}_hp2_straglr.tsv")
    params: 
        out_dir = os.path.join(config["OUT_DIR"], "straglr", "{sid}_{seqtype}_straglr")
        ref_fasta = cong["REF_FASTA"]
    threads: 1
    resources: mem_mb=2000
    log: "./logs/straglr.log"
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
        snp = os.path.join(config["OUT_DIR"], "1kG_SNP", "1kGP_high_coverage_Illumina.{chr}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"),
        tr_snp = os.path.join(config["OUT_DIR"], "TR_SNP_panel", "{chr}_TR_SNP_merged.vcf.gz")
        common_samples = rule.make_sample_info_file.output.common_info
    output: os.path.join(config["OUT_DIR"], "imputed_TR_SNPs", "{sid}_{seqtype}")
    params:
        out_dir = os.path.join(config["OUT_DIR"], "straglr", "{sid}_{seqtype}_straglr")
        beagle_jar = configs["BEAGLE_JAR_PATH"]
    threads: 2
    resources: mem_mb=4000
    log: "./logs/impute.log"
    shell:
        """
        mkdir -p {params.out_dir}/imputed_TRs
        java -Xmx4g -jar {params.beagle_jar} \
            gt=$(bcftools view -S {input.common_samples} {input.snp})  \
            ref=$(bcftools view -S {input.common_samples} {input.tr_snp}) \
            out={output}
        """
        
