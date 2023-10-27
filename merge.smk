configfile: "config.yaml"
from py.helpers import *

sample_dict = get_sample_dict(config, init=False)
sample_dict_expanded = get_sample_dict(config, init=True)

rule all:
    input:
        expand(os.path.join(config["DATA_DIR"], "hifi", "{sample}.hifi{ext}"),
            zip,
            sample=get_samples(sample_dict, "hifi"),
            num=get_num(sample_dict, "hifi"),
            ext=get_ext(sample_dict, "hifi")
        )

rule merge_short:
    output: os.path.join(config["DATA_DIR"], "short", "{sample}.short{ext}")
    params:
        sample=lambda wildcards: wildcards.sample,
        ext=lambda wildcards: wildcards.ext,
        files2merge=list_to_string(list(sample_dict_expanded["short"]["url"].keys())),
    conda: "envs/sam.yaml"
    shell: 
        '''
        # shouldnt need to merge short just rename it
        mv {params.files2merge} {output}
        '''

rule merge_hifi:
    output: os.path.join(config["DATA_DIR"], "hifi", "{sample}.hifi{ext}")
    params:
        sample=lambda wildcards: wildcards.sample,
        ext=lambda wildcards: wildcards.ext,
        files2merge=list_to_string(list(sample_dict_expanded["hifi"]["url"].keys())),
    conda: "envs/sam.yaml"
    shell: 
        '''
        samtools merge -o {output} {params.files2merge}
        '''

rule merge_assembly:
    output: os.path.join(config["DATA_DIR"], "assembly", "{sample}.assembly{ext}")
    params:
        sample=lambda wildcards: wildcards.sample,
        ext=lambda wildcards: wildcards.ext,
        files2merge=list_to_string(list(sample_dict_expanded["assembly"]["url"].keys())),
    conda: "envs/sam.yaml"
    shell: 
        '''
        # shouldnt need to merge assembly just rename it
        mv {params.files2merge} {output}
        '''