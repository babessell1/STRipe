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

# Define the number of cores and memory for each Snakemake run
declare -A resources
resources["download.smk"]="1,1000"
resources["index.smk"]="1,4000"
resources["call_trgt.smk"]="1,32000"

# Function to run Snakemake
run_snakemake() {
    snakefile="$1"
    cores_memory=(${resources[$snakefile]})
    cores="${cores_memory[0]}"
    mem_mb="${cores_memory[1]}"
    
    snakemake -s "$snakefile" -c "$config_file" --cores "$cores" --resources "mem_mb=$mem_mb"
}

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
    ext=$(get_ext "$sample_dict" "hifi")
    if [ "$datatype" != "HIFI" ]; then
        continue
    fi

    # Create a temporary manifest for the current sample
    temp_manifest="temp_manifest.csv"
    echo "$sample,$haplotype,$file_num,$datatype,$url" > "$temp_manifest"

    # Run each Snakefile with the specified configuration
    for snakefile in "${snakefiles[@]}"; do
        run_snakemake "$snakefile"
    done

    # Delete the temporary manifest and HIFI BAM file
    cleanup "$temp_manifest" "$sample" ".bam"
done < "$manifest_file"


