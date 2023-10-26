from helpers import *
configfile: "config.yaml"

sample_dict = get_sample_dict(configfile)

rule all:
    input:
        expand("data/short/{sample}.{datatype}", sample=wildcards.short, datatype=sample_dict["short"][wildcards.short]["datatype"])


rule download_short:
    # output should be in config["raw_dir""]
    output:
        "{raw_dir}/short_reads/{sample}.{datatype}}"
    params:
        raw_dir=config["raw_dir"],
        sample=lambda wildcards: wildcards.sample,
        datatype=lambda wildcards: sample_dict["short"][wildcards.sample]["datatype"],
        num=lambda wildcards: sample_dict["short"][wildcards.sample]["file_num"],
        url=lambda wildcards: sample_dict["short"][wildcards.sample]["url"]
    shell:
    
        """
        # if url is not s3 use wget
        mkdir -p raw_data/short_reads
        if [[ {params.url} != https://s3* ]]; then
            wget -O {output} {params.url}
        else
            aws s3 cp {params.url} {output}
        fi
        """

rule download_hifi:
    output:
        "{raw_dir}/hifi/{sample}.{datatype}"
    params:
        raw_dir=config["raw_dir"],
        sample=lambda wildcards: wildcards.sample,
        datatype=lambda wildcards: sample_dict["hifi"][wildcards.sample]["datatype"],
        url=lambda wildcards: sample_dict["hifi"][wildcards.sample]["url"]
    shell:
        """
        mkdir -p raw_data/hifi
        mkdir -p raw_data/assemblies
        # if url is not s3 use wget
        if [[ {params.url} != https://s3* ]]; then
            wget -O {output} {params.url}
        else
            aws s3 cp {params.url} {output}
        fi
        """
    
