import os 

configfile: "config.yaml"
from py.helpers import *

sample_dict = get_sample_dict(config, init=False)

rule_all = []

if config["PROCESS_SHORT"]:
    rule_all.extend([
        expand(os.path.join(config["DATA_DIR"], "short", "{sample}.short{ext}.{iext}"),
            zip,
            sample=get_samples(sample_dict, "short"),
            num=get_num(sample_dict, "short"),
            ext=get_ext(sample_dict, "short"),
            iext=get_iext(sample_dict, "short")
        )
    ])

if config["PROCESS_HIFI"]:
    rule_all.extend([
        expand(os.path.join(config["DATA_DIR"], "hifi", "{sample}.hifi{ext}.{iext}"),
            zip,
            sample=get_samples(sample_dict, "hifi"),
            num=get_num(sample_dict, "hifi"),
            ext=get_ext(sample_dict, "hifi"),
            iext=get_iext(sample_dict, "hifi")
        )
    ])

if config["PROCESS_ASSEMBLY"]:
    rule_all.extend([
        expand(os.path.join(config["DATA_DIR"], "assemblies", "{sample}.assembly{ext}.{iext}"),
            zip,
            sample=get_samples(sample_dict, "assembly"),
            num=get_num(sample_dict, "assembly"),
            ext=get_ext(sample_dict, "assembly"),
            iext=get_iext(sample_dict, "assembly")
        )
    ])


rule all:
    input: rule_all


rule get_short_index:
    input:  os.path.join(config["DATA_DIR"], "short", "{sample}.short{ext}")
    output: os.path.join(config["DATA_DIR"], "short", "{sample}.short{ext}.{iext}")
    params:
        url=lambda wildcards: sample_dict["short"]["url"][wildcards.sample],
        index=lambda wildcards: sample_dict["short"]["iext"][wildcards.sample]
    resources:
        mem_mb=2000
    threads: 1
    shell:
        '''
        wget -O "{output}" "{params.url}.{params.index}" || samtools index "{input}"
        '''


rule get_hifi_index:
    input: os.path.join(config["DATA_DIR"], "hifi", "{sample}.hifi{ext}")
    output: os.path.join(config["DATA_DIR"], "hifi", "{sample}.hifi{ext}.{iext}")
    params:
        url=lambda wildcards: sample_dict["hifi"]["url"][wildcards.sample],
        index=lambda wildcards: sample_dict["hifi"]["iext"][wildcards.sample]
    resources:
        mem_mb=4000
    threads: 1
    shell:
        '''
        # if url is not S3 use wget
        mkdir -p raw_data/short_reads
        if [[ ! "{params.url}" == "https://s3"* ]]; then
            wget -O "{output}" "{params.url}.{params.index}" || samtools index "{input}"
        else
            # Convert the URL to S3 format and download using AWS CLI
            s3_key=$(echo "{params.url}" | sed -e 's~https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=~~')
            aws s3 cp "s3://human-pangenomics/${{s3_key}}.{params.index}" "{output}" || samtools index "{input}"
        fi
        '''


rule get_assembly_index:
    input: os.path.join(config["DATA_DIR"], "assemblies", "{sample}.assembly{ext}")
    output: os.path.join(config["DATA_DIR"], "assemblies", "{sample}.assembly{ext}.{iext}")
    params:
        url=lambda wildcards: sample_dict["assemblies"]["url"][wildcards.sample],
        index=lambda wildcards: sample_dict["assemblies"]["iext"][wildcards.sample]
    resources:
        mem_mb=2000
    threads: 1
    shell:
        '''
        # if url is not S3 use wget
        mkdir -p raw_data/assemblies
        samtools faidx "{input}"
        '''