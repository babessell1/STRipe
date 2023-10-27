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
        wget -O "{output}" {params.url}
        '''


rule download_hifi:
    output: os.path.join(config["DATA_DIR"], "hifi", "{sample}.hifi{ext}")
    params:
        num=lambda wildcards: sample_dict["hifi"]["file_num"][wildcards.sample],
        url=lambda wildcards: sample_dict["hifi"]["url"][wildcards.sample]
    shell:
        '''
        # if url is not S3 use wget
        mkdir -p raw_data/short_reads
        if [[ ! "{params.url}" == "https://s3"* ]]; then
            wget -O "{output}" "{params.url}"
        else
            # Convert the URL to S3 format and download using AWS CLI
            s3_key=$(echo "{params.url}" | sed -e 's~https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=~~')
            aws s3 cp "s3://human-pangenomics/${{s3_key}}" "{output}"
        fi
        '''


rule get_short_index:
    input: rules.download_short.output
    output: os.path.join(config["DATA_DIR"], "short", "{sample}.short{ext}", "{iext}")
    params:
        num=lambda wildcards: sample_dict["short"]["file_num"][wildcards.sample],
        url=lambda wildcards: sample_dict["short"]["url"][wildcards.sample],
        ext=lambda wildcards: sample_dict["short"]["ext"][wildcards.sample],
        iext=lambda wildcards: sample_dict["short"]["iext"][wildcards.sample]
    conda: "envs/sam.yaml"
    shell:
        '''
        # if url is not S3 use wget
        mkdir -p raw_data/short_reads
        # try to download reads from same folder as cram/bam if it exists (just add index)
        wget -O "{output}" "{params.url}.{params.iext}" || samtools index "{input}"  
        '''

rule get_hifi_index:
    input: rules.download_hifi.output
    output: os.path.join(config["DATA_DIR"], "hifi", "{sample}.hifi{ext}", "{iext}")
    params:
        num=lambda wildcards: sample_dict["hifi"]["file_num"][wildcards.sample],
        url=lambda wildcards: sample_dict["hifi"]["url"][wildcards.sample],
        ext=lambda wildcards: sample_dict["hifi"]["ext"][wildcards.sample],
        iext=lambda wildcards: sample_dict["hifi"]["iext"][wildcards.sample]
    conda: "envs/sam.yaml"
    shell:
        '''
        # if url is not S3 use wget
        mkdir -p raw_data/short_reads
        # try to download reads from same folder as cram/bam if it exists (just add index)
        if [[ ! "{params.url}" == "https://s3"* ]]; then
            wget -O "{output}" "{params.url}.{params.iext}" || samtools index "{input}"
        else
            # Convert the URL to S3 format and download using AWS CLI
            s3_key=$(echo "{params.url}" | sed -e 's~https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=~~')
            aws s3 cp "s3://human-pangenomics/${{s3_key}}.{params.iext}" "{output}" || samtools index "{input}"
        fi
        '''