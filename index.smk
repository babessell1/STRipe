configfile: "config.yaml"
from py.helpers import *

sample_dict = get_sample_dict(config, init=True)

rule all:
    input:
        expand(os.path.join(config["DATA_DIR"], "short", "{sample}.short{ext}.{iext}"),
            zip,
            sample=get_samples(sample_dict, "short"),
            num=get_num(sample_dict, "short"),
            ext=get_ext(sample_dict, "short"),
            iext=get_iext(sample_dict, "short")
        ),
        expand(os.path.join(config["DATA_DIR"], "hifi", "{sample}.hifi{ext}.{iext}"),
            zip,
            sample=get_samples(sample_dict, "hifi"),
            num=get_num(sample_dict, "hifi"),
            ext=get_ext(sample_dict, "hifi"),
            iext=get_iext(sample_dict, "hifi")
        )


rule get_short_index:
    input:  os.path.join(config["DATA_DIR"], "short", "{sample}.short{ext}")
    output: os.path.join(config["DATA_DIR"], "short", "{sample}.short{ext}.{iext}")
    params:
        url=lambda wildcards: sample_dict["short"]["url"][wildcards.sample],
        index=lambda wildcards: sample_dict["short"]["iext"][wildcards.sample]
    conda: "envs/sam.yaml"
    shell:
        '''
        sample = "{sample}"
        # if url is not s3 use wget
        mkdir -p raw_data/short_reads
        wget -O "{output}" "{params.url}.{params.index}" || samtools index "{input}"
        '''


rule get_hifi_index:
    input: os.path.join(config["DATA_DIR"], "hifi", "{sample}.short{ext}")
    output: os.path.join(config["DATA_DIR"], "hifi", "{sample}.hifi{ext}.{iext}")
    params:
        url=lambda wildcards: sample_dict["hifi"]["url"][wildcards.sample],
        index=lambda wildcards: sample_dict["hifi"]["iext"][wildcards.sample]
    conda: "envs/sam.yaml"
    shell:
        '''
        # if url is not S3 use wget
        mkdir -p raw_data/short_reads
        # try to download reads from same folder as cram/bam if it exists (just add index)
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
    conda: "envs/sam.yaml"
    shell:
        '''
        # if url is not S3 use wget
        mkdir -p raw_data/assemblies
        wget -O "{output}" "{params.url}.{params.index}" || samtools faidx "{input}"
        '''
