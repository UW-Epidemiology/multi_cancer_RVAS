# Setup #

library(pacman)

p_load(tidyverse, data.table, openxlsx)

# Load and prep UKB #
UKB_single_results = NULL
for (i in c("ectoderm","mesoderm","endoderm","hormone","infectious","smoking")) {
  single_results = NULL
  for (j in 1:22) {
    chr_results = fread(paste0("/lindstroem/austin_working/Dissertation/UKBB/gene_based_results/cancer_groups/sparse_GRM/chr",j,"_",i,"_SAIGE_results.txt.singleAssoc.txt")) %>% 
      filter((AF_Allele2 < 0.01 | AF_Allele2 > 0.99) == TRUE, Allele2 != "UR") %>%
      mutate(group = i)
    if(j != 22 & nrow(chr_results)<1){
      next
    }
    single_results = rbind.data.frame(single_results,chr_results)
  }
  if(nrow(single_results)<1){
    next
  }
  UKB_single_results = rbind.data.frame(UKB_single_results,single_results)
}

## Check rare variant counts by analysis ##
UKB_single_results %>%
  group_by(group) %>%
  summarise(N = n())

## Apply P-value filter for significant results ##
UKB_single_results = UKB_single_results %>%
  filter(p.value < 5E-8)

annotations = NULL
for (i in unique(UKB_single_results$CHR) %>% as.integer()) {
  chr_annotations = fread(paste0("/lindstroem/austin_working/Dissertation/UKBB/variant_sets/filtered_gene_set_variants_to_QC/filtered_chr",i,".tsv.gz")) %>%
    select(c(Uploaded_variation,Gene,SYMBOL,Consequence)) %>%
    filter(Uploaded_variation %in% UKB_single_results$MarkerID)
  annotations = rbind.data.frame(annotations,chr_annotations)
}

UKB_single_results = merge(UKB_single_results,
                           annotations,
                           by.x = "MarkerID",
                           by.y="Uploaded_variation")

# Significant UKB burden test result tables 

## File output path
out_file <- "./austin_working/Dissertation/Meta/Aim1/scripts/result_summaries/tables/UKB_sig_single_tables.xlsx"

  table <- UKB_single_results %>%
    mutate(
      OR_CI = paste0(
        round(exp(BETA), 2), " (",
        round(exp(BETA+(SE*qnorm(0.025))), 2), ", ",
        round(exp(BETA+(SE*qnorm(0.975))), 2), ")"
      )
    ) %>%
    ungroup() %>%
    mutate(P = formatC(p.value, format = "E", digits = 2)) %>%
    arrange(SYMBOL, group) %>%
    select(MarkerID, SYMBOL, group, OR_CI, P, Consequence, AC_Allele2)
  
  write.xlsx(table,out_file)


  # Load and prep AoU #
  convert_variant_format <- function(variant) {
    variant <- gsub("^chr", "", variant)  # Remove "chr" from the beginning
    variant <- gsub(":", "-", variant)    # Replace ":" with "-"
    return(variant)
  }
  
  AoU_single_results = NULL
  for (i in c("ectoderm","mesoderm","endoderm","hormone","infectious","smoking")) {
    single_results = NULL
    for (j in 1:22) {
      chr_results = fread(paste0("/lindstroem/austin_working/Dissertation/AoU/gene_based_results/cancer_groups/sparse_GRM/chr",j,"_",i,"_SAIGE_results.txt.singleAssoc.txt")) %>%
        filter((AF_Allele2 < 0.01 | AF_Allele2 > 0.99) == TRUE, Allele2 != "UR") %>%
        mutate(group = i)
      if(j != 22 & nrow(chr_results)<1){
        next
      }
      single_results = rbind.data.frame(single_results,chr_results)
    }
    if(nrow(single_results)<1){
      next
    }
    AoU_single_results = rbind.data.frame(AoU_single_results,single_results)
  }
  
## Check rare variant counts by analysis ##
  AoU_single_results %>%
    group_by(group) %>%
    summarise(N = n())
  
## Apply P-value filter for significant results ##
  AoU_single_results = AoU_single_results %>%
    filter(p.value < 5E-8)  
  
## Convert marker format for annotation merge ##  
  AoU_single_results = AoU_single_results %>%
    mutate(Uploaded_variation = convert_variant_format(MarkerID))
  
  annotations = NULL
  for (i in unique(AoU_single_results$CHR) %>% as.integer()) {
    chr_annotations = fread(paste0("/lindstroem/austin_working/Dissertation/AoU/variant_sets/filtered_gene_set_variants_to_QC/filtered_chr",i,".tsv.gz")) %>%
      select(c(Uploaded_variation,Gene,SYMBOL,Consequence)) %>%
      filter(Uploaded_variation %in% AoU_single_results$Uploaded_variation)
    annotations = rbind.data.frame(annotations,chr_annotations)
  }
  
  AoU_single_results = merge(AoU_single_results, annotations,
                             by="Uploaded_variation")
  
  # Significant AoU burden test result tables
  
  ## File output path
  out_file <- "./austin_working/Dissertation/Meta/Aim1/scripts/result_summaries/tables/AoU_sig_single_tables.xlsx"
  
  table <- AoU_single_results %>%
    mutate(
      OR_CI = paste0(
        round(exp(BETA), 2), " (",
        round(exp(BETA+(SE*qnorm(0.025))), 2), ", ",
        round(exp(BETA+(SE*qnorm(0.975))), 2), ")"
      )
    ) %>%
    ungroup() %>%
    mutate(P = formatC(p.value, format = "E", digits = 2)) %>%
    arrange(SYMBOL, group) %>%
    select(MarkerID, SYMBOL, group, OR_CI, P, Consequence, AC_Allele2)
  
  write.xlsx(table,out_file)
  
  
  