# Set up #
rm(list = ls())
library(tidyverse)
library(data.table)
library(stringr)
source("/lindstroem/UKBB_AoU_WGS_annotations/scripts/VEP_annotation_utilites.R")
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

# Loop over chromosomes 1 through 22
for (chr in 1:22) {
  
  # Load data #
  input_file <- paste0("/lindstroem/UKBB_AoU_WGS_annotations/UKBB/WES_annotations/WES_variants_chr", chr, "_VEP_ANN.tsv.gz")
  data <- read_VEP(input_file)
  
  # Apply filters #
  
  ## 1: Only variants annotated to MANE Select transcripts ##
  chr_data <- data %>% filter(is.na(MANE_SELECT) == FALSE)
  
  ## 2: Only coding and splice site variants ##
  chr_data <- chr_data %>% filter(is.na(chr_data$EXON) == FALSE |
                                    sapply(Consequence, function(x) any(str_detect(x, splice_site_variants))))
  
  ## 3A. Variants with CADD >20 (excluding synonymous variants) ##
  chr_high_CADD_data <- chr_data %>% filter(CADD_PHRED > 20 &
                                              strsplit(as.character(Consequence), ",") %!in% "synonymous_variant")
  
  ## 3B. Synonymous variants ##
  chr_synonymous_data <- chr_data %>% filter(strsplit(as.character(Consequence), ",") %in% "synonymous_variant")
  
  ## 3C. Variants with missing CADD scores ##
  chr_missing_CADD_data <- chr_data %>% filter(is.na(CADD_PHRED) == TRUE)
  
  chr_missing_CADD_data_include <- chr_missing_CADD_data %>%
    filter(sapply(Consequence, function(x) any(str_detect(x, pLoF_variants))) |
             sapply(CLIN_SIG, function(x) any(str_detect(x, clin_pathogenic))) |
             sapply(am_class, function(x) any(str_detect(x, "likely_pathogenic"))) |
             REVEL > 0.5 |
             sapply(SIFT, function(x) any(str_detect(x, SIFT_deleterious))) |
             sapply(PolyPhen, function(x) any(str_detect(x, Poly_Phen_deleterious))))
  
  # Create final filtered dataset #
  filtered_data <- rbind(chr_high_CADD_data, chr_synonymous_data, chr_missing_CADD_data_include)
  filtered_data <- filtered_data %>% arrange(Location)
  
  # Write final filtered dataset #
  output_file <- paste0("/lindstroem/austin_working/Dissertation/UKBB/filtered_gene_set_variants_to_QC/filtered_chr", chr, ".tsv.gz")
  fwrite(filtered_data, output_file, compress = "gzip", row.names = FALSE, quote = FALSE, sep = "\t")
  
  # Cleanup #
  rm(list = c("chr_high_CADD_data", "chr_missing_CADD_data_include", "chr_missing_CADD_data", "chr_data", "chr_synonymous_data"))
  
  cat("Chromosome", chr, "processed successfully.\n")
}