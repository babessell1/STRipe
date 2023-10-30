configfile: "config.yaml"
from py.helpers import *

sample_dict = get_sample_dict(config, init=False)

rule all:
    input:

rule get_pileup:
    input: os.path.join(config["DATA_DIR"], "short", "{sample}.short{ext}")
    output: ??????
    params:
        ref=config["REF_FASTA"]
        prefix=os.path.join(config["DATA_DIR"], "ehdn", "{sample}")
    conda: "envs/sam.yaml"
    shell: 
        '''
        ExpansionHunter --reads {input} \
                --reference {params.ref} \
                --variant-catalog {params.catalog_json} \
                --output-prefix {params.prefix} \
                --analysis-mode streaming
        '''