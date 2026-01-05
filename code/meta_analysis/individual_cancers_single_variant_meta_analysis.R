# Setup #
library(pacman)

p_load(rlang, tidyverse, data.table,
       R.utils, metafor, meta)

# Define useful functions #
convert_annotation_id_to_AoU_format <- function(variant) {
  # Add "chr" if not already present
  variant <- ifelse(grepl("^chr", variant), variant, paste0("chr", variant))
  # Replace "-" with ":"
  variant <- gsub("-", ":", variant)
  return(variant)
}

## Identify variants in common ##
find_common_variants_dt <- function(df1, df2, window_size = 0) {
  # Convert to data.table
  dt1 <- as.data.table(df1)
  dt2 <- as.data.table(df2)
  
  # Ensure correct types
  dt1[, CHR := as.integer(as.character(CHR))]
  dt2[, CHR := as.integer(as.character(CHR))]
  dt1[, POS := as.numeric(as.character(POS))]
  dt2[, POS := as.numeric(as.character(POS))]
  
  # Add global row index
  dt1[, global_idx1 := .I]
  dt2[, global_idx2 := .I]
  
  # Define function to get MAJ/MIN pair
  get_major_minor <- function(a1, a2, af2) {
    maj <- ifelse(af2 > 0.5, a2, a1)
    min <- ifelse(af2 > 0.5, a1, a2)
    paste(maj, min, sep = "/")
  }
  
  dt1[, AllelePair := get_major_minor(Allele1, Allele2, AF_Allele2)]
  dt2[, AllelePair := get_major_minor(Allele1, Allele2, AF_Allele2)]
  
  # Initialize result indices
  matched_idx1 <- integer()
  matched_idx2 <- integer()
  common_ids <- character()
  
  common_chrs <- intersect(unique(dt1$CHR), unique(dt2$CHR))
  
  for (chr in common_chrs) {
    chr_dt1 <- dt1[CHR == chr]
    chr_dt2 <- dt2[CHR == chr]
    
    if (window_size == 0) {
      # Exact match
      m <- merge(
        chr_dt1[, .(CHR, POS, AllelePair, global_idx1)],
        chr_dt2[, .(CHR, POS, AllelePair, global_idx2)],
        by = c("CHR", "POS", "AllelePair")
      )
    } else {
      # Windowed overlap
      chr_dt1[, start := POS - window_size]
      chr_dt1[, end := POS + window_size]
      chr_dt2[, POS2 := POS]
      
      setkey(chr_dt1, CHR, start, end)
      
      m <- foverlaps(
        chr_dt2[, .(CHR, POS, POS2, AllelePair, global_idx2)],
        chr_dt1[, .(CHR, start, end, POS, AllelePair, global_idx1)],
        by.x = c("CHR", "POS", "POS"),
        by.y = c("CHR", "start", "end"),
        nomatch = 0
      )[AllelePair == i.AllelePair]
    }
    
    if (nrow(m) > 0) {
      matched_idx1 <- c(matched_idx1, m$global_idx1)
      matched_idx2 <- c(matched_idx2, m$global_idx2)
      common_ids <- c(common_ids, paste0("chr", m$CHR, "_", m$POS, "_", m$AllelePair))
    }
  }
  
  # Extract filtered rows
  dt1_common <- dt1[matched_idx1]
  dt2_common <- dt2[matched_idx2]
  
  dt1_common[, Common_ID := common_ids]
  dt2_common[, Common_ID := common_ids]
  
  return(list(df1_common = dt1_common, df2_common = dt2_common))
}


## Harmonize to minor allele effects ##
harmonize_effects_to_minor_allele <- function(common_variants) {
  dt1_all <- as.data.table(common_variants$df1_common)
  dt2_all <- as.data.table(common_variants$df2_common)
  
  # Ensure CHR is treated as character
  dt1_all[, CHR := as.character(CHR)]
  dt2_all[, CHR := as.character(CHR)]
  
  chromosomes <- intersect(unique(dt1_all$CHR), unique(dt2_all$CHR))
  results_list <- vector("list", length(chromosomes))
  
  for (i in seq_along(chromosomes)) {
    chr <- chromosomes[i]
    
    # Copy subsets to avoid modifying originals
    dt1 <- copy(dt1_all[CHR == chr])
    dt2 <- copy(dt2_all[CHR == chr])
    
    keep_vars <- c("SYMBOL", "CHR", "POS", "MarkerID", "Allele1", "Allele2",
                   "AC_Allele2", "AF_Allele2", "Common_ID", "Consequence",
                   "BETA", "SE", "p.value", "N_case","N_ctrl")
    
    dt1 <- dt1[, ..keep_vars]
    dt2 <- dt2[, ..keep_vars]
    
    # Rename columns
    setnames(dt1, old = keep_vars[!(keep_vars %in% c("CHR", "POS", "Common_ID"))],
             new = paste0(keep_vars[!(keep_vars %in% c("CHR", "POS", "Common_ID"))], "_UKB"))
    setnames(dt2, old = keep_vars[!(keep_vars %in% c("CHR", "POS", "Common_ID"))],
             new = paste0(keep_vars[!(keep_vars %in% c("CHR", "POS", "Common_ID"))], "_AoU"))
    
    merged <- merge(dt1, dt2, by = c("CHR", "POS", "Common_ID"))
    
    merged[, Flip_UKB := AF_Allele2_UKB > 0.5]
    merged[, Flip_AoU := AF_Allele2_AoU > 0.5]
    
    merged[, BETA_UKB := ifelse(Flip_UKB, -BETA_UKB, BETA_UKB)]
    merged[, BETA_AoU := ifelse(Flip_AoU, -BETA_AoU, BETA_AoU)]
    
    merged[, MAF_UKB := pmin(AF_Allele2_UKB, 1 - AF_Allele2_UKB)]
    merged[, MAF_AoU := pmin(AF_Allele2_AoU, 1 - AF_Allele2_AoU)]
    
    results_list[[i]] <- merged
  }
  
  final_result <- rbindlist(results_list)
  return(final_result)
}

# Load variant annotations #
UKB_annotations = NULL
for (i in 1:22) {
  chr_variants = fread(paste0("/lindstroem/austin_working/Dissertation/UKBB/variant_sets/final_variant_sets/chr",i,"_final_variants.txt")) %>%
    unlist() %>%
    as.vector()
  chr_annotations = fread(paste0("/lindstroem/austin_working/Dissertation/UKBB/variant_sets/filtered_gene_set_variants_to_QC/filtered_chr",i,".tsv.gz")) %>%
    select(c(Uploaded_variation,Gene,SYMBOL,Consequence)) %>%
    filter(Uploaded_variation %in% chr_variants)
  UKB_annotations = rbind.data.frame(UKB_annotations,chr_annotations)
}

AoU_annotations = NULL
for (i in 1:22) {
  chr_annotations = fread(paste0("/lindstroem/austin_working/Dissertation/AoU/variant_sets/filtered_gene_set_variants_to_QC/filtered_chr",i,".tsv.gz")) %>%
    select(c(Uploaded_variation,Gene,SYMBOL,Consequence))
  AoU_annotations = rbind.data.frame(AoU_annotations,chr_annotations)
}

AoU_annotations = AoU_annotations %>%
  mutate(Uploaded_variation = convert_annotation_id_to_AoU_format(Uploaded_variation))


# Run meta-analyses by cancer type #

meta_genes = c("ATM", "CHEK2", "BRCA2", "BRCA1", "CDKN2A", "TP53", "HAL", "HOXB13", "MSH6",  
               "OCA2", "DNMT3A", "MLH1", "PALB2", "CPVL", "FLG", "MITF", "PPM1D", "PRDM7", "RTEL1",  
               "SRSF2", "ASXL1", "FAM111A", "IFIH1", "NF1", "POT1", "TET2", "PLEKHA4", "HLA-DPB1",  
               "MICA", "JAK2", "IDH2", "IGLL5","SAMHD1")

results = NULL
  
for (i in c("anus", "bladder", "bone", "breast",
            "brain", "cervix", "colorectal",
            "endometrium", "esophagus", "HL",
            "kidney", "leukemias", "liver",
            "lung", "melanoma", "myeloma",
            "NHL", "neck", "oral",
            "ovary", "pancreas", "prostate",
            "stomach", "testis", "thyroid",
            "eye","non_melanoma")) {
  

  
  UKBB = fread(paste0("/lindstroem/austin_working/Dissertation/UKBB/gene_based_results/individual_cancers_combined/",i,"_SAIGE_results.txt.singleAssoc.txt")) %>%
    filter((N_case_het + N_case_hom) > 0)
  
  AoU = fread(paste0("/lindstroem/austin_working/Dissertation/AoU/gene_based_results/individual_cancers_combined/",i,"_SAIGE_results.txt.singleAssoc.txt")) %>%
    filter((N_case_het + N_case_hom) > 0)
  
  
  UKB_results = merge(UKBB,
                      UKB_annotations,
                      by.x = "MarkerID",
                      by.y = "Uploaded_variation",
                      all.x = TRUE) %>%
    filter(SYMBOL %in% meta_genes)
  
  AoU_results = merge(AoU,
                      AoU_annotations,
                      by.x = "MarkerID",
                      by.y = "Uploaded_variation",
                      all.x = TRUE) %>%
    filter(SYMBOL %in% meta_genes)
  
  ## Harmonize variants ##
  common = find_common_variants_dt(UKB_results,AoU_results)
  
  harmonized <- harmonize_effects_to_minor_allele(common) %>%
    mutate(CHR = as.integer(CHR)) %>%
    arrange(CHR, POS)
  
  harmonized_max_MAF = harmonized %>%
    filter(MAF_UKB <= 0.01 & MAF_AoU <= 0.01)
  
  ## Perform meta-analysis ##
  meta_results <- harmonized_max_MAF %>%
    transmute(
      Gene = SYMBOL_UKB,
      CHR,
      POS,
      ID = Common_ID,
      Consequence = Consequence_UKB,
      MAF_UKB,
      MAF_AoU,
      N_case_UKB,
      N_case_AoU,
      N_control_UKB = N_ctrl_UKB,
      N_control_AoU = N_ctrl_AoU,
      UKB_BETA = BETA_UKB,
      AoU_BETA = BETA_AoU,
      UKB_SE = SE_UKB,
      AoU_SE = SE_AoU,
      P_UKB = p.value_UKB,
      P_AoU = p.value_AoU
    )
  
  meta_stats <- pmap_dfr(
    list(harmonized_max_MAF$BETA_UKB, harmonized_max_MAF$BETA_AoU,
         harmonized_max_MAF$SE_UKB, harmonized_max_MAF$SE_AoU),
    function(b1, b2, se1, se2) {
      fit <- rma(yi = c(b1, b2), sei = c(se1, se2), method = "FE")
      tibble(
        BETA = as.numeric(fit$beta),
        BETA_SE = fit$se,
        BETA_CI_LB = fit$ci.lb,
        BETA_CI_UB = fit$ci.ub,
        P = fit$pval,
        I2 = fit$I2
      )
    }
  )
  
  meta_results <- bind_cols(meta_results, meta_stats) %>%
    mutate(Cancer = i)
  
  results = rbind.data.frame(results,meta_results)
}


fwrite(results, 
       paste0("/lindstroem/austin_working/Dissertation/Meta/Aim1/individual_cancers/individual_cancers_single_variant_meta.txt"),
       compress = "gzip")
