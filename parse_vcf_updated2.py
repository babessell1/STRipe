import os
import pandas as pd
import json
import gzip


def extract_variant_info(json_data):
    reference_regions = []
    confidence_intervals = []
    flanking_reads_list = []
    inrepeat_reads_list = []
    spanning_reads_list = []
    
    locus_results = json_data.get('LocusResults', {})

    for locus_id, locus_info in locus_results.items():
        variants = locus_info.get('Variants', {})

        #extract the information we need from JSON file
        for variant_id, variant_info in variants.items():
            reference_region = variant_info.get('ReferenceRegion')
            confidence_interval = variant_info.get('GenotypeConfidenceInterval')
            flanking_reads = variant_info.get('CountsOfFlankingReads')
            inrepeat_reads = variant_info.get('CountsOfInrepeatReads')
            spanning_reads = variant_info.get('CountsOfSpanningReads')

            if reference_region:
                modified_chr = reference_region.replace(':', '_').replace('-', '_')
                reference_regions.append(modified_chr)
                confidence_intervals.append(confidence_interval)
                flanking_reads_list.append(flanking_reads)
                inrepeat_reads_list.append(inrepeat_reads)
                spanning_reads_list.append(spanning_reads)
                
    return reference_regions, confidence_intervals, flanking_reads_list, inrepeat_reads_list, spanning_reads_list

def parse_vcf_and_json(sample_id, vcf_directory, json_directory):
    vcf_file = f"{sample_id}.hifi.sorted.vcf.gz"
    json_file = f"{sample_id}.json"

    vcf_file_path = os.path.join(vcf_directory, vcf_file)
    json_file_path = os.path.join(json_directory, json_file)

    if os.path.exists(vcf_file_path) and os.path.exists(json_file_path):
        # parse JSON file
        with open(json_file_path, 'r') as json_file:
            json_data = json.load(json_file)
        reference_regions, confidence_intervals, flanking_reads_list, inrepeat_reads_list, spanning_reads_list = extract_variant_info(json_data)
        
       # parse VCF file
        chrom_data = []
        with gzip.open(vcf_file_path, "rt") as file:
            for line in file:
                if line.startswith("#CHROM"):
                    continue
                parts = line.strip().split("\t")
                if len(parts) >= 8:
                    chrom = parts[0]
                    pos = parts[1]
                    info = parts[7].split(";")
                    trid = info[0].split("=")[1]
                    end_pos = info[1].split("=")[1]
                    motif = info[2].split("=")[1]
                    len_motif = len(motif)
                    struc = info[3].split("=")[1]
                    MC = parts[9].split(":")[4]
                    # Check if trid is in the reference_regions
                    if trid in reference_regions:
                    ## Find the index of trid in the reference_regions list
                        idx = reference_regions.index(trid)
                        
                        # Append additional information to chrom_data list
                        chrom_data.append([
                            chrom, pos, trid, end_pos, motif, len_motif, struc, MC,
                            confidence_intervals[idx],
                            flanking_reads_list[idx],
                            inrepeat_reads_list[idx],
                            spanning_reads_list[idx]
                        ])

        # create a df
        column_names = ['chrom', 'pos', 'trid', 'end_pos', 'motif', 'len_motif', 'struc', 'MC', 
                    'confidence_interval', 'flanking_reads', 'inrepeat_reads', 'spanning_reads']
        
        vcf_df = pd.DataFrame(chrom_data, columns=column_names)
        # Extract sample ID from the file name and use it as the sample name
        sample_id_from_filename = sample_id.split(".")[0]
        vcf_df['sample_name'] = sample_id_from_filename
        
        region_length_list = []

        # calculating the repeat region length for every row and adding it as an additional column in df
        for row in vcf_df["trid"]:
            first_pos = int(row.split("_")[1])
            sec_pos = int(row.split("_")[2])
            region_len = sec_pos - first_pos
            region_length_list.append(region_len)

        vcf_df["len_repeat_region"] = region_length_list  # add length of repeat region as an extra column

        # filter for repeat regions <= 10 and allocate that to a new df
        vcf_data_10 = vcf_df[(vcf_df["len_repeat_region"] >= 10)]

        return vcf_df, vcf_data_10
    else:
        print(f"Error: Either {vcf_file} or {json_file} is missing.")

####################################################################################################
#applying above functions to vcf files in directory for trgt output
####################################################################################################
vcf_directory = "/nfs/turbo/dcmb-class/bioinf593/groups/group_05/output/trgt"
json_directory = "/nfs/turbo/dcmb-class/bioinf593/groups/group_05/raw/eh"

# Get a list of all VCF files in the directory
vcf_files = [f.split('.')[0] for f in os.listdir(vcf_directory) if f.endswith('.vcf.gz')]

# Initialize an empty DataFrame to store the combined results
combined_df = pd.DataFrame()
combined_10_df = pd.DataFrame()

# Loop through each VCF file
for sample_id in vcf_files:
    vcf_df, vcf_10_df = parse_vcf_and_json(sample_id, vcf_directory, json_directory)
    if vcf_df is not None:
        combined_df = combined_df.append(vcf_df, ignore_index=True)
        combined_10_df = combined_10_df.append(vcf_10_df, ignore_index=True)

# Export the combined DataFrame
output_csv_path = '/home/vivianra/all_vcf_parsed.csv'
combined_df.to_csv(output_csv_path, index=False, sep="\t")
print(f"Combined DataFrame exported to {output_csv_path}")

output_csv_10_path = '/home/vivianra/all_vcf_10_combined.csv'
combined_10_df.to_csv(output_csv_10_path, index=False, sep="\t")
print(f"Combined vcf_10_df exported to {output_csv_10_path}")

