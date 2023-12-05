# Parallel Pileup Parse
# Auden Bahr


# IF THE ORDER OF THE COLUMNS chr, pos and end_pos in the vcf, or any in the pileup CHANGES, THE FUNCTIONS NEED TO BE EDITED

require(data.table)
require(parallel)

path_to_pileups <- '/nfs/turbo/dcmb-class/bioinf593/groups/group_05/pileup_temp/'
sample_list <- list.files(path = path_to_pileups, pattern = "\\.cut.pileup$") # 41 - this is the input list

# start of operations on each sample
# define for now

parse_sample <- function(sample_file) {

# run once with either file
#path_to_vcf <- '/nfs/turbo/dcmb-class/bioinf593/groups/group_05/output/trgt/repeatregion_10_parsedvcf.txt'
path_to_vcf <- '/nfs/turbo/dcmb-class/bioinf593/groups/group_05/output/trgt/repeatregion_all_parsedvcf.txt'
path_to_pileups <- '/nfs/turbo/dcmb-class/bioinf593/groups/group_05/pileup_temp/'

output_path <- '/nfs/turbo/dcmb-class/bioinf593/groups/group_05/output/depth/all_repeat_regions'

sample <- gsub("\\.cut.pileup$", "", sample_file)

read_vcf_chr <- sprintf(
  "grep --text '%s' %s | cut -f1 | uniq",
  sample,
  paste(path_to_vcf)
)

sample_chr <- data.table::fread(cmd = read_vcf_chr, header = FALSE)
sample_chr_vec <- sample_chr$V1

for (chromosome in sample_chr_vec) { # change column order
  print(chromosome)
  
  if (!file.exists(paste0(output_path, sample, '_', chromosome, '_depth.txt'))) {
  # command to read in the position column 
  # filter for chromosome
  read_pileup_pos <- sprintf(
    "grep --text '%s' %s | awk '{print $2}'",
    chromosome,
    paste0(path_to_pileups, sample_file)
  )
  
  # command to read in start position of repeat from vcf
  # filter for sample and chromosome
  read_vcf_start <- sprintf(
    "grep --text '%s' %s | grep --text '%s' | awk '{print $2}'",
    sample,
    paste(path_to_vcf),
    chromosome
  )
  
  # command to read in end position of repeat from vcf
  # filter for sample and chromosome
  read_vcf_end <- sprintf(
    "grep --text '%s' %s | grep --text '%s' | awk '{print $4}'",
    sample,
    paste(path_to_vcf),
    chromosome
  )
  
  read_repeat_id <- sprintf(
    "grep --text '%s' %s | grep --text '%s' | awk '{print $3}'",
    sample,
    paste('/nfs/turbo/dcmb-class/bioinf593/groups/group_05/output/trgt/repeatregion_10_parsedvcf.txt'),
    chromosome
  )
  
  # command to read in pileup depth
  # filter for chromosome
  read_pileup_depth <- sprintf(
    "grep --text '%s' %s | awk '{print $4}'",
    chromosome,
    paste0(path_to_pileups, sample_file)
  )
  
  # load positions in pileup file
  pileup_pos <- data.table::fread(cmd = read_pileup_pos, header = FALSE)
  
  # load start positions in vcf file
  vcf_start <- data.table::fread(cmd = read_vcf_start, header = FALSE)
  
  setkey(pileup_pos, V1) 
  setkey(vcf_start, V1)
  
  # get row index of start position in pileup file
  start_index <- pileup_pos[vcf_start, which = TRUE, on = .(V1 = V1)]
  
  # remove vcf start positions
  rm(vcf_start)
  
  # load repeat end positions from vcf
  vcf_end <- data.table::fread(cmd = read_vcf_end, header = FALSE)
  setkey(vcf_end, V1)
  
  # get row index of end position in pileup file
  end_index <- pileup_pos[vcf_end, which = TRUE, on = .(V1 = V1)]
  
  # remove vcf end positions and pileup positions
  rm(vcf_end)
  rm(pileup_pos)
  
  # read in repeat id
  trid <- data.table::fread(cmd = read_repeat_id, header = FALSE)
  
  # check if there are any NAs in start_index or end_index
  if (any(is.na(start_index)) || any(is.na(end_index))) {
    # get the indices of the NAs
    na_indices_start <- which(is.na(start_index))
    na_indices_end <- which(is.na(end_index))
    
    # combine the indices
    all_na_indices <- unique(c(na_indices_start, na_indices_end))
    print(all_na_indices)
    # remove the NA values in start_index, end_index, and repeat_id
    start_index <- start_index[-all_na_indices]
    end_index <- end_index[-all_na_indices]
    trid <- trid[-all_na_indices, ]
  }
  
  # read in pileup depth
  pileup_depth <- data.table::fread(cmd = read_pileup_depth)
  
  # range: 250 upstream and downstream (500 values total)
  range=250-1
  
  # slice pileup depth
  result <- mapply(function(start, end)
    
    #             [row indices, column index=1]
    c(pileup_depth[(max(1, start-range)):start, 1],
      pileup_depth[end:(min(end+range, nrow(pileup_depth))), 1]),
    start_index, end_index, SIMPLIFY=FALSE
    )
  
  # flatten result into dataframe with comma-separated depth values
  df_flat <- data.frame(
    depth_values = sapply(result, function(x) paste(unlist(x), collapse = ', '))
  )
  
  # remove result of slicing and pileup depth values
  rm(result)
  rm(pileup_depth)
  
  # create output file with sample, repeat id, and depth values
  trid_w_depth <- cbind(sample, trid, df_flat)
  
  # write output file
  write.table(
    trid_w_depth,
    paste0(output_path, sample, '_', chromosome, '_depth.txt'),
    row.names=FALSE, col.names = FALSE, quote = FALSE)
  
  print(sample)
  
  }
}

}

cores <- detectCores()
r <- mclapply(sample_list, parse_sample, mc.cores = cores)


