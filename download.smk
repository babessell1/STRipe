from py.helpers import *
configfile: "config.yaml"

sample_dict = get_sample_dict(config)

print(sample_dict)

rule all:
    input:
        expand(os.path.join(config["DATA_DIR"], "short", "{sample}.short{ext}"),
            sample=list(sample_dict["short"]["url"].keys()),
            num=[val for val in list(sample_dict["short"]["file_num"].values())],
            ext=[val for val in list(sample_dict["short"]["ext"].values())]
        ),
        expand(os.path.join(config["DATA_DIR"], "hifi", "{sample}.hifi{ext}"),
            sample=list(sample_dict["hifi"]["url"].keys()),
            num=[val for val in list(sample_dict["hifi"]["file_num"].values())],
            ext=[val for val in list(sample_dict["hifi"]["ext"].values())]
        )


rule download_short:
    output: os.path.join(config["DATA_DIR"], "short", "{sample}.short{ext}")
    params:
        num=lambda wildcards: sample_dict["short"]["file_num"][wildcards.sample],
        url=lambda wildcards: sample_dict["short"]["url"][wildcards.sample]
    shell:
        '''
        # if url is not s3 use wget
        mkdir -p raw_data/short_reads
        if [[ ! {params.url} == "https://s3"* ]]; then
            wget -O {output} {params.url}
        else
            aws s3 cp {params.url} {output}
        fi
        '''

rule download_hifi:
    output: os.path.join(config["DATA_DIR"], "hifi", "{sample}.hifi{ext}")
    params:
        num=lambda wildcards: sample_dict["hifi"]["file_num"][wildcards.sample],
        url=lambda wildcards: sample_dict["hifi"]["url"][wildcards.sample]
    shell:
        '''
        mkdir -p raw_data/hifi
        mkdir -p raw_data/assemblies
        # if url is not s3 use wget
        if [[ ! {params.url} == "https://s3"* ]]; then
            wget -O {output} {params.url}
        else
            aws s3 cp {params.url} {output}
        fi
        '''
    
