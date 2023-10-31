configfile: "config.yaml"
from py.helpers import *

sample_dict = get_sample_dict(config, init=False)

rule all:
    input:
        expand(os.path.join(config["OUT_DIR"], "trgt", "{sample}.?????????????????"),
            zip,
            sample=get_samples(sample_dict, "hifi")
        )

rule call_trgt:
    input: os.path.join(config["DATA_DIR"], "hifi", "{sample}.hifi{ext}")
    output: os.path.join(config["OUT_DIR"], "trgt", "{sample}.?????????????????")
    params:
        ref=config["REF_FASTA"],
        trgt=config["TRGT"],
        trgt_bed=config["TRGT_BED"],
    resources:
        mem_mb=32000
    threads: 16
    conda: "envs/trgt.yaml"
    shell:
        '''
        
        '''

