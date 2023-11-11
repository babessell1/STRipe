#!/bin/bash

#SBATCH --account=bioinf593f23_class 
#SBATCH --job-name=run_trgt
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --mem=32GB
#SBATCH --ntasks=1
#SBATCH --time=7:59:00
#SBATCH --output=logs/trgt.out
#SBATCH --error=logs/trgt.err


# run one by one to prevent disk space issues

MODE=$1

# set cores to ntasks
CORES="$SLURM_NTASKS"

manifest_file="manifests/hifi_manifest.csv"

if [ "$MODE" == "inline" ]; then
    # Read the manifest file and process each line
    while IFS=',' read -r sample haplotype file_num datatype url; do

        if [ "$datatype" != "HIFI" ]; then
            continue
        fi

        # Create a temporary manifest for the current sample
        temp_manifest="./manifests/temp_manifest.csv"
        echo "sample_name,haplotype,file_num,datatype,long_read_url" > "$temp_manifest"
        echo "$sample,$haplotype,$file_num,$datatype,$url" >> "$temp_manifest"
        hifi_bam_file="/nfs/turbo/dcmb-class/bioinf593/groups/group_05/raw/hifi/${sample}.hifi.bam"
        trgt_file="/nfs/turbo/dcmb-class/bioinf593/groups/group_05/output/trgt/${sample}.hifi.sorted.vcf.gz"
        realign_bam_file="/nfs/turbo/dcmb-class/bioinf593/groups/group_05/output/trgt/${sample}.hifi.sorted.spanning.bam"

        # if trgt file does not exist, then run the pipeline
        if [ ! -f "$trgt_file" ]; then
            # Run each Snakefile with the specified configuration
            snakemake -s "download.smk" -c "$config_file" --cores 1 --resources "mem_mb=1000"
            snakemake -s "index.smk" -c "$config_file" --cores 1 --resources "mem_mb=4000"
            snakemake -s "call_trgt.smk" -c "$config_file" --cores 1 --resources "mem_mb=32000"

            # Delete hifi bam file and temp_manifest and spanning relaligned bam
            rm "$hifi_bam_file"
            rm "$realign_bam_file"
        fi

        rm "$temp_manifest"
        
    done < "$manifest_file"

elif [ "$MODE" == "download" ]; then
    mem=$(($CORES * 1000))
    snakemake -s "download.smk" --cores "$CORES" --resources "mem_mb=${mem}" -n -r
elif [ "$MODE" == "index" ]; then
    mem=$(($CORES * 4000))
    snakemake -s "index.smk" --cores "$CORES" --resources "mem_mb=${mem}" -n -r
elif [ "$MODE" == "call_trgt" ]; then
    mem=$(($CORES * 32000))
    snakemake -s "call_trgt.smk" --cores "$CORES" --resources "mem_mb=${mem}" -n -r
fi

