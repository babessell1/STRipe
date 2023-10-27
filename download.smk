from py.helpers import *
configfile: "config.yaml"

sample_dict = get_sample_dict(config)

rule all:
    input:
        expand("data/short/{sample}.short.{ext}", sample=wildcards.short, datatype=sample_dict["short"]["ext"][wildcards.short]),
        expand("data/hifi/{sample}.hifi.{ext}", sample=wildcards.hifi, datatype=sample_dict["hifi"]["ext"][wildcards.hifi])


rule download_short:
    # output should be in config["raw_dir""]
    output:
        "{raw_dir}/short/{sample}.short.{ext}"
    params:
        raw_dir=config["raw_dir"],
        sample=lambda wildcards: wildcards.sample,
        ext=lambda wildcards: sample_dict["short"]["ext"][wildcards.sample],
        num=lambda wildcards: sample_dict["short"]["file_num"][wildcards.sample],
        url=lambda wildcards: sample_dict["short"]["url"][wildcards.sample]
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
        datatype=lambda wildcards: sample_dict["hifi"]["ext"][wildcards.sample],
        url=lambda wildcards: sample_dict["hifi"]["url"][wildcards.sample]
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
    
