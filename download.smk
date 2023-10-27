from py.helpers import *
configfile: "config.yaml"

sample_dict = get_sample_dict(config)

rule all:
    input:
        expand("data/short/{sample}.short.{ext}", sample=wildcards.short, datatype=sample_dict["short"]["ext"][wildcards.short]),
        expand("data/hifi/{sample}.hifi.{ext}", sample=wildcards.hifi, datatype=sample_dict["hifi"]["ext"][wildcards.hifi])


rule download_short:
    # output should be in config["raw_dir""]
    input: lambda wildcards: "{raw_dir}/touch/{sample}.short.touch"
    output: "{raw_dir}/short/{sample}.short.{ext}"
    params:
        root_dir=config["ROOT_DIR"],
        raw_dir=config["RAW_DIR"],
        sample=lambda wildcards: wildcards.sample,
        ext=lambda wildcards: sample_dict["short"]["ext"][wildcards.hifi],
        num=lambda wildcards: sample_dict["short"]["file_num"][wildcards.hifi],
        url=lambda wildcards: sample_dict["short"]["url"][wildcards.short]
    shell:
        """
        # if url is not s3 use wget
        mkdir -p raw_data/short_reads
        if [[ {input.url} != https://s3* ]]; then
            wget -O {output} {params.url}
        else
            aws s3 cp {input.url} {output}
        fi
        """

rule download_hifi:
    input: lambda wildcards: "{raw_dir}/touch/{sample}.hifi.touch"
    output: "{raw_dir}/hifi/{sample}.{datatype}"
    params:
        root_dir=config["ROOT_DIR"],
        raw_dir=config["RAW_DIR"],
        sample=lambda wildcards: wildcards.sample,
        datatype=lambda wildcards: sample_dict["hifi"]["ext"][wildcards.sample],
        url=lambda wildcards: sample_dict["short"]["url"][wildcards.sample]
    shell:
        """
        mkdir -p raw_data/hifi
        mkdir -p raw_data/assemblies
        # if url is not s3 use wget
        if [[ {input.url} != https://s3* ]]; then
            wget -O {output} {params.url}
        else
            aws s3 cp {input.url} {output}
        fi
        """
    
