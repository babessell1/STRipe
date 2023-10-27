from snakemake import snakemake

# Import the other Snakefiles
include: "download.smk"
include: "index.smk"
include: "merge.smk"
include: "call_trgt.smk"

# Define a rule to run all the other rules
rule run_all:
    input:
        rules.download_hifi.output,
        rules.get_hifi_index.output
    output:
        # Define a dummy output to track the completion of all rules
        touch("all_rules_completed.txt")
    run:
        with open("manifests/long_manifest.csv") as handle:
            for line in handle:
                sample, haplotype, file_num, datatype, url = line.split(",")
                if datatype != "HIFI":
                    continue
                # create a temp file to run each smk on its own mini manifest
                with open("temp_manifest.csv", "w") as temp_handle:
                    temp_handle.write(line)

                # run the snakefile on the temp manifest
                snakemake(
                    snakefile="download.smk",
                    configfile="config.yaml",
                    cores=1,
                )
                snakemake(
                    snakefile="index.smk",
                    configfile="config.yaml",
                    cores=1,
                )
                snakemake(
                    snakefile="merge.smk",
                    configfile="config.yaml",
                    cores=1,
                )
                snakemake(
                    snakefile="call_trgt.smk",
                    configfile="config.yaml",
                    cores=1,
                )
                

