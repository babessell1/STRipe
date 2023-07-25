import os
import sys
from warnings import warn
from py.helpers import *
configfile: "config.yaml"

print(get_seqtype(os.path.join(config["OUT_DIR"], "sample_info_common.tsv")))
'''
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
        "tricolor",  
        "{sid}_{seqtype}_tricolor.chk"
    ), zip,
    sid=get_samp_id(os.path.join(config["OUT_DIR"], "sample_info_common.tsv")),
    seqtype=get_seqtype(os.path.join(config["OUT_DIR"], "sample_info_common.tsv"))
),
'''
rule_all = [
    #expand(
    #    os.path.join(config["OUT_DIR"],
    #        "straglr", 
    #        "{sid}_{seqtype}_straglr",
    #        "{sid}_{seqtype}_hp1_straglr.tsv"
    #    ), zip,
    #    sid=get_samp_id(os.path.join(config["OUT_DIR"], "sample_info_common.tsv")),
    #    seqtype=get_seqtype(os.path.join(config["OUT_DIR"], "sample_info_common.tsv"))
    #),
    #expand(
    #    os.path.join(config["OUT_DIR"],
    #        "straglr", 
    #        "{sid}_{seqtype}_straglr",
    #        "{sid}_{seqtype}_hp2_straglr.tsv"
    #    ), zip,
    #    sid=get_samp_id(os.path.join(config["OUT_DIR"], "sample_info_common.tsv")),
    #    seqtype=get_seqtype(os.path.join(config["OUT_DIR"], "sample_info_common.tsv"))
    #),
    #expand(
    #    os.path.join(config["OUT_DIR"],
    #        "imputed_TR_SNPs",
    #        "{chr}"
    #    ),
    #    chr=get_chromosomes(config["CHROMOSOMES"]),
    #)
    expand(
        os.path.join(config["OUT_DIR"],
            "vamos", 
            "{sid}_{seqtype}_vamos",
            "{sid}_{seqtype}_hp1.vcf"
        ), zip,
        sid=get_samp_id(os.path.join(config["OUT_DIR"], "sample_info_common.tsv")),
        seqtype=get_seqtype(os.path.join(config["OUT_DIR"], "sample_info_common.tsv"))
    ),
    expand(
        os.path.join(config["OUT_DIR"],
            "vamos", 
            "{sid}_{seqtype}_vamos",
            "{sid}_{seqtype}_hp2.vcf"
        ), zip,
        sid=get_samp_id(os.path.join(config["OUT_DIR"], "sample_info_common.tsv")),
        seqtype=get_seqtype(os.path.join(config["OUT_DIR"], "sample_info_common.tsv"))
    )

]


rule all: input: rule_all

rule vamos:
    input:
        hp1 = os.path.join(config["DATA_DIR"], "{sid}.{seqtype}.h1.sorted.bam"),
        hp2 = os.path.join(config["DATA_DIR"], "{sid}.{seqtype}.h2.sorted.bam")
    output:
        hp1 = os.path.join(config["OUT_DIR"], "vamos", "{sid}_{seqtype}_vamos", "{sid}_{seqtype}_hp1.vcf"),
        hp2 = os.path.join(config["OUT_DIR"], "vamos", "{sid}_{seqtype}_vamos", "{sid}_{seqtype}_hp2.vcf")
    
    params:
        sid = "{sid}",
        seqtype = "{seqtype}",
        emotifs = config["EMOTIFS"],
        vamos = config["VAMOS_PATH"]
    threads: 8
    resources: mem_mb = 30000
    conda: "envs/vamos.yml"
    shell:
        """
        {params.vamos} --contig -b {input.hp1} -r {params.emotifs} -s {params.sid} -o {output.hp1} -t {threads}
        {params.vamos} --contig -b {input.hp2} -r {params.emotifs} -s {params.sid} -o {output.hp2} -t {threads}
        """



'''
rule tricolor:
    output: 
        sensor = directory(os.path.join(config["OUT_DIR"], "tricolor", "sensor", "{sid}_{seqtype}_sensor")),
        refer = directory(os.path.join(config["OUT_DIR"], "tricolor", "refer", "{sid}_{seqtype}_refer")),
        check = os.path.join(config["OUT_DIR"], "tricolor", "{sid}_{seqtype}_tricolor.chk")
    input:
        hp1 = os.path.join(config["DATA_DIR"], "{sid}.{seqtype}.h1.sorted.bam"),
        hp2 = os.path.join(config["DATA_DIR"], "{sid}.{seqtype}.h2.sorted.bam")
    params: 
        data_dir = os.path.join(config["OUT_DIR"]),
        ref_fasta = config["REF_FASTA"],
        sid = "{sid}",
        seqtype = "{seqtype}"
    threads: 4
    resources: mem_mb=30000
    log: "logs/tricolor_{sid}_{seqtype}.log"
    shell:
        """
        . ../mambaforge/etc/profile.d/conda.sh
        #conda init bash
        conda activate tricolorenv
        mkdir -p {params.data_dir}/tricolor/sensor/{params.sid}_{params.seqtype}_sensor
        mkdir -p {params.data_dir}/tricolor/refer

        TRiCoLOR SENSoR -bam {input.hp1} {input.hp2} -o {output.sensor} >> {log}
        TRiCoLOR REFER -g {params.ref_fasta} -bam {input.hp1} {input.hp2} -bed {output.sensor}/*.bed -o {output.refer} >> {log}
        conda deactivate
        touch {output.check}
        """
'''

#rule straglr:
#    input:
#        hp1 = os.path.join(config["DATA_DIR"], "{sid}.{seqtype}.h1.sorted.bam"),
#        hp2 = os.path.join(config["DATA_DIR"], "{sid}.{seqtype}.h2.sorted.bam")
#    output: 
#        hp1 = os.path.join(config["OUT_DIR"], "straglr", "{sid}_{seqtype}_straglr", "{sid}_{seqtype}_hp1_straglr.tsv"),
#        hp2 = os.path.join(config["OUT_DIR"], "straglr", "{sid}_{seqtype}_straglr", "{sid}_{seqtype}_hp2_straglr.tsv")
#    params: 
#        out_dir = os.path.join(config["OUT_DIR"], "straglr", "{sid}_{seqtype}_straglr"),
#        ref_fasta = config["REF_FASTA"],
#        straglr_path = config["STRAGLR_PATH"]
#    threads: 32
#    resources: mem_mb=15000
#    log: "logs/straglr.{sid}.{seqtype}.log"
#    shell:
#        """
#        . ../mambaforge/etc/profile.d/conda.sh
#        #conda init bash
#        conda activate straglr
#        mkdir -p {params.out_dir}/strglr_temp
#        python {params.straglr_path} \
#            {input.hp1} \
#            {params.ref_fasta} \
#            "$(cut -f 1 -d "." <<< "{output.hp1}")" \
#            --tmpdir {params.out_dir}/strglr_temp \
#            --chroms chr21
#        >> {log}
#        
#        python {params.straglr_path} \
#            {input.hp2} \
#            {params.ref_fasta} \
#            "$(cut -f 1 -d "." <<< "{output.hp2}")"  \
#            --tmpdir {params.out_dir}/strglr_temp \
#            --chroms chr21
#        >> {log}
#        conda deactivate
#        """


#rule impute:
#    input: 
#        snp = os.path.join(config["SNP_DIR"], "1kGP_high_coverage_Illumina.{chr}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"),
#        tr_snp = os.path.join(config["SNP_TR_PANEL_DIR"], "{chr}_TR_SNP_merged.vcf.gz"),
#        common_samples = os.path.join(config["OUT_DIR"], "sample_info_common.tsv")
#    output: os.path.join(config["OUT_DIR"], "imputed_TR_SNPs", "{chr}")
#    params:
#        out_dir = os.path.join(config["OUT_DIR"]),
#        beagle_jar = config["BEAGLE_JAR_PATH"],
#        chr = "{chr}"
#    threads: 16
#    resources: mem_mb=122000
#    log: "logs/impute.{chr}.log"
#    shell:
#        """
#        . ../mambaforge/etc/profile.d/conda.sh
#        #conda init bash
#        conda activate stromboli
#        mkdir -p {params.out_dir}/imputed_TR_SNPs
#        mkdir -p {params.out_dir}/impute_tmp
#        cat {input.common_samples} | sed 's/|/ /' | awk '{{print $1}}' | uniq > {params.out_dir}/impute_tmp/{params.chr}.tmp 
#        bcftools view -S {params.out_dir}/impute_tmp/{params.chr}.tmp {input.tr_snp} -Oz -o {params.out_dir}/impute_tmp/{params.chr}_gt.vcf.gz
#        java -Xmx122g -jar {params.beagle_jar} \
#            gt={params.out_dir}/impute_tmp/{params.chr}_gt.vcf.gz \
#            ref={input.tr_snp} \
#            out={output} \
#        >> {log}
#        rm {params.out_dir}/impute_tmp/{params.chr}.tmp
#        conda deactivate
#        """