import os

configfile: "config.yaml"
from py.helpers import *

sample_dict = get_sample_dict(config, init=False)

rule all:
    input:
        expand(os.path.join(config["OUT_DIR"], "trgt", "{sample}.hifi.sorted.vcf.gz"),
            zip,
            sample=get_samples(sample_dict, "hifi")
        )

rule call_trgt:
    input: os.path.join(config["DATA_DIR"], "hifi", "{sample}.hifi.bam")
    output: os.path.join(config["OUT_DIR"], "trgt", "{sample}.hifi.sorted.vcf.gz")
    params:
        ref=config["REF_FASTA"],
        trgt=config["TRGT_PATH"],
        trgt_bed=config["TRGT_BED"],
    resources:
        mem_mb=32000
    threads: 1
    shell:
        '''
        # prefix is output without .tar.gz
        prefix=$(echo {output} | sed 's/.vcf.gz//g')
        {params.trgt} --genome {params.ref} \
            --repeats {params.trgt_bed} \
            --reads {input} \
            --output-prefix $prefix \
            --threads {threads}
        '''

    
