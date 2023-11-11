#!/bin/bash

# take pileup file name as input
pileup_file=$1

# take range as input
range=$2

declare -A vcf_info
tail -n +2 /nfs/turbo/dcmb-class/bioinf593/groups/group_05/pileup_temp/vcf_parsed_HG00438.txt | while read -r chr start_pos end_pos
do
    vcf_info["$chr,$start_pos,$end_pos"]=1
done

awk -v range="$range" -v OFS=',' '{
    for (key in vcf_info) {
        split(key, vcf_data, ",")
        if ($1 == vcf_data[1] && $2 >= vcf_data[2] - range && $2 <= vcf_data[3] + range) {
            print $1, $2, $4
        }
    }
}' "$pileup_file" > depth_col.csv

paste -d' ' sample_vcf.txt depth_col.csv > combined_sample.csv
