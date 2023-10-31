configfile: "config.yaml"
from py.helpers import *

sample_dict = get_sample_dict(config, init=False)
sample_dict_expanded = get_sample_dict(config, init=True)

rule_all = []

if config["PROCESS_SHORT"] == True:
    rule_all.extend([
        expand(os.path.join(config["DATA_DIR"], "short", "{sample}.short{ext}"),
            zip,
            sample=get_samples(sample_dict, "short"),
            num=get_num(sample_dict, "short"),
            ext=get_ext(sample_dict, "short")
        )
    ])
if config["PROCESS_HIFI"] == True:
    rule_all.extend([
        expand(os.path.join(config["DATA_DIR"], "hifi", "{sample}.hifi{ext}"),
            zip,
            sample=get_samples(sample_dict, "hifi"),
            num=get_num(sample_dict, "hifi"),
            ext=get_ext(sample_dict, "hifi")
        )
    ])
if config["PROCESS_ASSEMBLY"] == True:
    rule_all.extend([
        expand(os.path.join(config["DATA_DIR"], "assembly", "{sample}.assembly{ext}"),
            zip,
            sample=get_samples(sample_dict, "assembly"),
            num=get_num(sample_dict, "assembly"),
            ext=get_ext(sample_dict, "assembly")
        )
    ])


rule all:
    input: rule_all


if config["PROCESS_SHORT"] == True:
    rule merge_short:
        output: os.path.join(config["DATA_DIR"], "short", "{sample}.short{ext}")
        params:
            sample=lambda wildcards: wildcards.sample,
            ext=lambda wildcards: wildcards.ext,
            files2merge=list_to_string(list(sample_dict_expanded["short"]["url"].keys())),
        conda: "envs/sam.yaml"
        resources:  # 1GB
            mem_mb=1000
        threads: 1
        shell: 
            '''
            # shouldnt need to merge short just rename it
            mv {params.files2merge} {output}
            # remove input
            rm {params.files2merge}
            '''


if config["PROCESS_HIFI"] == True:
    rule merge_hifi:
        output: os.path.join(config["DATA_DIR"], "hifi", "{sample}.hifi{ext}")
        params:
            sample=lambda wildcards: wildcards.sample,
            ext=lambda wildcards: wildcards.ext,
            files2merge=list_to_string(list(sample_dict_expanded["hifi"]["url"].keys())),
        conda: "envs/sam.yaml"
        resources:
            mem_mb=1000
        threads: 1
        shell: 
            '''
            samtools merge -o {output} {params.files2merge}
            rm {params.files2merge}
            # sort
            #samtools sort -o {output} {output}
            '''


if config["PROCESS_ASSEMBLY"] == True:
    rule merge_assembly:
        output: os.path.join(config["DATA_DIR"], "assembly", "{sample}.assembly{ext}")
        params:
            sample=lambda wildcards: wildcards.sample,
            ext=lambda wildcards: wildcards.ext,
            files2merge=list_to_string(list(sample_dict_expanded["assembly"]["url"].keys())),
        conda: "envs/sam.yaml"
        resources:
            mem_mb=1000
        threads: 1
        shell: 
            '''
            # shouldnt need to merge assembly just rename it
            mv {params.files2merge} {output}
            rm {params.files2merge}
            '''