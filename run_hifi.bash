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

# Define the configuration file
config_file="config.yaml"

# Function to delete the temporary manifest and HIFI BAM file
cleanup() {
    temp_manifest="$1"
    sample="$2"
    ext="$3"
    
    rm "$temp_manifest"
    hifi_bam_file="$DATA_DIR/hifi/${sample}.hifi.bam"
    if [ -e "$hifi_bam_file" ]; then
        rm "$hifi_bam_file"
    fi
}

# Read the manifest file and process each line
while IFS=',' read -r sample haplotype file_num datatype url; do

    if [ "$datatype" != "HIFI" ]; then
        continue
    fi

    # Create a temporary manifest for the current sample
    temp_manifest="temp_manifest.csv"
    echo "$sample,$haplotype,$file_num,$datatype,$url" > "$temp_manifest"

    # Run each Snakefile with the specified configuration
    # Run each Snakefile with the specified configuration
    for snakefile in "${snakefiles[@]}"; do
        case "$snakefile" in
            "download.smk")
                snakemake -s "$snakefile" -c "$config_file" --cores 1 --resources "mem_mb=1000"
                ;;
            "index.smk")
                snakemake -s "$snakefile" -c "$config_file" --cores 1 --resources "mem_mb=4000"
                ;;
            "call_trgt.smk")
                snakemake -s "$snakefile" -c "$config_file" --cores 1 --resources "mem_mb=32000"
                ;;
        esac
    done

    # Delete the temporary manifest and HIFI BAM file
    cleanup "$temp_manifest" "$sample" ".bam"
done < "$manifest_file"


