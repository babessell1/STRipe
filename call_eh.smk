configfile: "config.yaml"
from py.helpers import *

sample_dict = get_sample_dict(config, init=False)

rule all:
    input:

rule call_eh:
    input: os.path.join(config["DATA_DIR"], "short", "{sample}.short{ext}")
    output:
        json=os.path.join(config["DATA_DIR"], "eh", "{sample}.json"),
        vcf=os.path.join(config["DATA_DIR"], "eh", "{sample}.vcf")
        re_bam=os.path.join(config["DATA_DIR"], "eh", "{sample}_realigned.bam"),
    params:
        ref=config["REF_FASTA"],
        prefix=os.path.join(config["DATA_DIR"], "eh", "{sample}"),
        eh=config["EH"]
    conda: "envs/sam.yaml"
    threads: 16
    shell:
        '''
        {params.eh} --reads {input} \
                --reference {params.ref} \
                --variant-catalog {params.catalog_json} \
                --output-prefix {params.prefix} \
                --analysis-mode streaming \
                --threads {threads}
        '''

rule subset_json:
    input: output.json
    output: os.path.join(config["DATA_DIR"], "eh", "{sample}_large.json")
    threads: 1
    resources:
        mem_mb=1000
    shell:
        '''
        jq '.LocusResults | with_entries(select(.value.Variants[].CountsOfInrepeatReads != "()"))' {input} > {output}
        '''
