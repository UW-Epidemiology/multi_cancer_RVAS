# Set up #
rm(list = ls())
library(tidyverse)
library(data.table)
library(stringr)
setDTthreads(percent = 100)

`%!in%` = Negate(`%in%`)

# Create some variant functional annotation category string variables #
splice_site_variants = c("splice_acceptor_variant", "splice_donor_variant")

pLoF_variants <- c("transcript_ablation", "splice_acceptor_variant",
                   "splice_donor_variant", "stop_gained", 
                   "frameshift_variant", "stop_lost",
                   "start_lost", "transcript_elongation",
                   "feature_elongation", "feature_truncation")

# Create some variant pathogenicity/delteriousness string variables #
clin_pathogenic = c("likely_pathogenic","pathogenic",
                    "pathogenic/likely_pathogenic",
                    "likely_pathogenic/pathogenic")

SIFT_deleterious = c("deleterious","deleterious_low_confidence")
Poly_Phen_deleterious = c("possibly_damaging","probably_damaging")

# Loop over chromosomes 1 to 22
for (chr in 1:22) {
  # Initialize filtered_data at the start of each chromosome
  filtered_data <- NULL
  part <- 1
  
  repeat {
    # Generate the file path dynamically
    file_path <- paste0("/projects/lindstroem/UKBB_AoU_WGS_annotations/AoU/v8_final_block_annotations/chr",chr,"/WGS_v8_variants_chr", chr, "_part", part, "_VEP_ANN.tsv.gz")
    
    # Check if the file exists to avoid errors
    if (!file.exists(file_path)) {
      if (part == 1) {
        cat("No data found for chromosome", chr, "\n")
      }
      break
    }
    
    # Load data
    data <- fread(file_path, na.strings = "-")
    
    # Step 1: Only variants annotated to MANE Select transcripts
    chr_data <- data %>% filter(!is.na(MANE_SELECT))
    
    # Skip part if chr_data is empty after filtering in Step 1
    if (nrow(chr_data) == 0) {
      cat("Chromosome", chr, "part", part, "has no data after Step 1 filtering and will be skipped.\n")
      part <- part + 1
      next
    }
    
    # Step 2: Only coding and splice site variants
    chr_data <- chr_data %>%
      mutate(splice_site_variant = sapply(Consequence, function(x) any(str_detect(x, splice_site_variants)))) %>%
      filter(!is.na(EXON) | splice_site_variant)
    
    # Skip part if chr_data is empty after filtering in Step 2
    if (nrow(chr_data) == 0) {
      cat("Chromosome", chr, "part", part, "has no data after Step 2 filtering and will be skipped.\n")
      part <- part + 1
      next
    }
    
    # Step 3A: Variants with CADD >20 (excluding synonymous variants)
    chr_high_CADD_data <- chr_data %>% filter(CADD_PHRED > 20 & !str_detect(as.character(Consequence), "synonymous_variant"))
    
    # Step 3B: Synonymous variants
    chr_synonymous_data <- chr_data %>% filter(str_detect(as.character(Consequence), "synonymous_variant"))
    
    # Step 3C: Variants with missing CADD scores
    chr_missing_CADD_data <- chr_data %>% filter(is.na(CADD_PHRED))
    
    if(nrow(chr_missing_CADD_data) > 0){
      chr_missing_CADD_data_include <- chr_missing_CADD_data %>%
        filter(sapply(Consequence, function(x) any(str_detect(x, pLoF_variants))) |
                 sapply(CLIN_SIG, function(x) any(str_detect(x, clin_pathogenic))) |
                 sapply(am_class, function(x) any(str_detect(x, "likely_pathogenic"))) |
                 REVEL > 0.5 |
                 sapply(SIFT, function(x) any(str_detect(x, SIFT_deleterious))) |
                 sapply(PolyPhen, function(x) any(str_detect(x, Poly_Phen_deleterious))))
    }else{
      chr_missing_CADD_data_include = NULL
    }
    

    # Combine filtered data from each part
    filtered_data_part <- rbind(chr_high_CADD_data, chr_synonymous_data, chr_missing_CADD_data_include)
    filtered_data_part <- filtered_data_part %>% arrange(Location)
    
    # Skip if no data in this part after combining filters
    if (nrow(filtered_data_part) == 0) {
      cat("Chromosome", chr, "part", part, "has no data after combining filters and will be skipped.\n")
      part <- part + 1
      next
    }
    
    # Append part data to chromosome-level data
    if (is.null(filtered_data)) {
      filtered_data <- filtered_data_part
    } else {
      filtered_data <- rbind(filtered_data, filtered_data_part)
    }
    
    # Cleanup for the current part
    rm(list = c("chr_high_CADD_data", "chr_synonymous_data", "chr_missing_CADD_data", "chr_missing_CADD_data_include", "chr_data", "data", "filtered_data_part"))
    
    cat("Chromosome", chr, "part", part, "processed successfully.\n")
    
    # Increment part number
    part <- part + 1
  }
  
  # Check if we have any data for this chromosome and write it
  if (!is.null(filtered_data) && nrow(filtered_data) > 0) {
    output_file_path <- paste0("/projects/lindstroem/austin_working/Dissertation/AoU/variant_sets/filtered_gene_set_variants_to_QC/filtered_chr", chr, ".tsv.gz")
    fwrite(filtered_data, output_file_path, compress = "gzip", row.names = FALSE, quote = FALSE, sep = "\t")
  } else {
    cat("Chromosome", chr, "has no data to save after processing all parts.\n")
  }
}