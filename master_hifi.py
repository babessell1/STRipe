#!/usr/bin/env python3

import snakemake
from snakemake.io import expand
from py.helpers import *
import os

# Define the Snakefiles you want to run
snakefiles = ["download.smk", "index.smk", "call_trgt.smk"]

# Define your manifest file
manifest_file = "manifests/hifi_manifest.csv"

# Define the configuration file
config_file = "config.yaml"

# Define the number of cores and memory for each Snakemake run
resources = {"download.smk": {"cores": 1, "mem_mb": 1000},
             "index.smk": {"cores": 1, "mem_mb": 4000},
             "call_trgt.smk": {"cores": 1, "mem_mb": 32000}}

sample_dict = get_sample_dict(config, init=False)

# Read the manifest file and process each line
with open(manifest_file) as handle:
    for line in handle:
        sample, haplotype, file_num, datatype, url = line.strip().split(",")
        ext=get_ext(sample_dict, "hifi")
        if datatype != "HIFI":
            continue

        # Create a temporary manifest for the current sample
        temp_manifest = f"temp_manifest.csv"
        with open(temp_manifest, "w") as temp_handle:
            temp_handle.write(line)

        # Run each Snakefile with the specified configuration
        for snakefile in snakefiles:
            snakemake(
                snakefile=snakefile,
                configfile=config_file,
                cores=resources[snakefile]["cores"],
                resources={"mem_mb": resources[snakefile]["mem_mb"]}
            )

        # Delete the temporary manifest
        os.remove(temp_manifest)

        # Delete the HIFI BAM file
        hifi_bam_file = os.path.join(config["DATA_DIR"], "hifi", f"{sample}.hifi{ext}")
        if os.path.exists(hifi_bam_file):
            os.remove(hifi_bam_file)
