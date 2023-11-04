from snakemake import snakemake as smk
configfile: "config.yaml"
from py.helpers import *
#import os

# Import the other Snakefiles
include: "download.smk"
include: "index.smk"
include: "merge.smk"
include: "call_trgt.smk"

sample_dict = get_sample_dict(config, init=False)

# Define a rule to run all the other rules
rule run_all:
    output:
        # Define a dummy output to track the completion of all rules
        "all_rules_completed.txt"
    run:
        with open("manifests/hifi_manifest.csv") as handle:
            for line in handle:
                sample, haplotype, file_num, datatype, url = line.split(",")
                ext = sample_dict["hifi"]["ext"][sample]
                if datatype != "HIFI":
                    continue
                # create a temp file to run each smk on its own mini manifest
                with open("temp_manifest.csv", "w") as temp_handle:
                    temp_handle.write(line)

                # run the snakefile on the temp manifest
                smk(
                    snakefile="download.smk",
                    configfile="config.yaml",
                    cores=1,
                    resources={"mem_mb": 1000}
                )
                smk(
                    snakefile="index.smk",
                    configfile="config.yaml",
                    cores=1,
                    resources={"mem_mb": 4000}
                )
                smk(
                    snakefile="call_trgt.smk",
                    configfile="config.yaml",
                    cores=1,
                    resources={"mem_mb": 32000}
                )
                # delete the temp manifest
                os.remove("temp_manifest.csv")
                # delete the hifi bam file
                os.remove(os.path.join(config["DATA_DIR"], "hifi", f"{sample}.hifi{ext}"))
                

