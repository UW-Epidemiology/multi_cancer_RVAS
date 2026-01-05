# Download exome filtered srWGS bim files (use Python kernel)#

# Setup (R kernel)#
rm(list = ls())
install.packages("pacman")
library(pacman)

p_load(tidyverse,data.table,R.utils)

# Create variant lists and group files #

## Define categories for annotation groups ##
pLoF_terms <- c("transcript_ablation", "splice_acceptor_variant", 
                "splice_donor_variant", "stop_gained", 
                "frameshift_variant", "stop_lost", 
                "start_lost", "transcript_elongation", 
                "feature_elongation", "feature_truncation")
missense_terms <- "missense_variant"
synonymous_terms <- "synonymous_variant"

## Run QC and group file creation loop ##
for (i in 1:22) {
  ## Load data ##
  chr_anno_vars = fread(paste0("./rare_variant_sets/initial_filtered_annotations/filtered_chr",i,".tsv.gz"),
                   select = c("Uploaded_variation","Gene","SYMBOL","Consequence"))
  chr_anno_vars$V2 <- paste0("chr",gsub("-", ":", chr_anno_vars$Uploaded_variation))
  chr_bim_vars = fread(paste0("./rare_variant_sets/PLINK_exome_bims/exome.chr",i,".bim"))
  chr_annotations = inner_join(chr_anno_vars,chr_bim_vars,by="V2") %>% select(c("V2","Gene","SYMBOL","Consequence"))
  
  ## Replace empty SYMBOL values with gene values ##
  chr_annotations <- chr_annotations %>%
    mutate(SYMBOL = if_else(is.na(SYMBOL) | SYMBOL == "", Gene, SYMBOL))
  
  ## Get distinct genes with at least 3 unique variants ##
  chr_genes <- chr_annotations %>%
    group_by(SYMBOL) %>%
    filter(n_distinct(V2) >= 3) %>% 
    ungroup() %>%
    distinct(SYMBOL) %>%
    pull() %>%
    na.omit()
  
  ## Filter for selected genes and ensure unique variants ##
  chr_annotations <- chr_annotations %>%
    filter(SYMBOL %in% chr_genes) %>%
    distinct(V2, .keep_all = TRUE)
  
  ## Variant annotations ##
  chr_annotations <- chr_annotations %>%
    mutate(var_ann = case_when(
      Consequence %in% pLoF_terms ~ "pLoF",
      Consequence %in% missense_terms ~ "missense",
      Consequence %in% synonymous_terms ~ "synonymous",
      TRUE ~ "other"
    ))
  
  ## Prepare output file ##
  output_file <- paste0("./rare_variant_sets/group_files/chr", i, "_group_file.txt")
  fileConn <- file(output_file, "wt")  
  
  for (j in chr_genes) {
    gene_variants <- chr_annotations %>%
      filter(SYMBOL == j)
    
    ## Check for duplicate variants within the gene ##
    duplicate_snps <- gene_variants %>%
      count(V2) %>%
      filter(n > 1)
    
    if (nrow(duplicate_snps) > 0) {
      warning(paste("Duplicate SNPs detected in gene", j, ":", paste(duplicate_snps$V2, collapse = ", ")))
    }
    
    ## Construct var and anno rows, ensuring no trailing spaces ##
    var_list <- paste(c(j, "var", gene_variants$V2), collapse = " ")
    ann_list <- paste(c(j, "anno", gene_variants$var_ann), collapse = " ")
    
    writeLines(trimws(var_list), fileConn)
    writeLines(trimws(ann_list), fileConn)
  }
  
  close(fileConn)  
  
  ## Write unique final variants per chromosome ##
  fwrite(chr_annotations %>% select(V2) %>% distinct(),
         paste0("./rare_variant_sets/final_variant_sets/chr", i, "_final_variants.txt"),
         col.names = FALSE, row.names = FALSE,
         quote = FALSE)
}

# Transfer files to the workspace bucket