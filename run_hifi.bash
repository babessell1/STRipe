#!/bin/bash

#SBATCH --account=bioinf593f23_class 
#SBATCH --job-name=stripe_download
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --mem=32GB
#SBATCH --ntasks=11
#SBATCH --time=6:00:00
#SBATCH --output=logs/hifi.out
#SBATCH --error=logs/hifi.err

#!/bin/bash

# Define the Snakefiles you want to run
snakefiles=("download.smk" "index.smk" "call_trgt.smk")

# Define your manifest file
manifest_file="manifests/hifi_manifest.csv"

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
    realign_bam_file="/nfs/turbo/dcmb-class/bioinf593/groups/group_05/output/trgt/${sample}.hifi.sorted.spanning.bam

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


