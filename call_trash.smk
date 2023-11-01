configfile: "config.yaml"
from py.helpers import *

sample_dict = get_sample_dict(config, init=False)

rule all:
    input:
        expand(os.path.join(config["DATA_DIR"], "trash", "{sample}.?????????????????"),
            zip,
            sample=get_samples(sample_dict, "assembly")
        )

rule call_trgt:
    input: os.path.join(config["DATA_DIR"], "short", "{sample}.assembly{ext}")
    output: os.path.join(config["OUT_DIR"], "trash", "{sample}.?????????????????")
    params:
        ref=config["REF_FASTA"],
        trash=config["TRASH"],
    shell:
        '''

        '''
