# Setup #
library(pacman)
p_load(data.table,tidyverse, flextable)

## Load data ##
single_results = fread("/lindstroem/austin_working/Dissertation/Meta/Aim1/individual_cancers/individual_cancers_single_variant_meta.txt") 



## Filter to nominally significant results ##
all_sig = single_results %>%
  filter(P < 0.05)

fwrite(all_sig,file = paste0(paste0("/lindstroem/austin_working/Dissertation/Meta/Aim1/individual_single_for_release/single_meta_individual_cancers_P_value_under_0.05.txt")))

## Filter to results with Bonferroni corrected significance #

bon_sig_tresh = 0.05 /length(unique(single_results$ID))

bon_all_sig = single_results %>%
  filter(P < bon_sig_tresh)
 
# Reformat for table #

formatted_sig_group_single_results = bon_all_sig %>%
  mutate(Meta_OR = paste0(round(exp(BETA),2),
                          " (",
                          round(exp(BETA_CI_LB),2),
                          " – ",
                          round(exp(BETA_CI_UB),2),
                          ")"),
         UKB_OR = paste0(round(exp(UKB_BETA),2),
                         " (",
                         round(exp(UKB_BETA + UKB_SE*qnorm(0.025)),2),
                         " – ",
                         round(exp(UKB_BETA + UKB_SE*qnorm(0.975)),2),
                         ")"),
         AoU_OR = paste0(round(exp(AoU_BETA),2),
                         " (",
                         round(exp(AoU_BETA + AoU_SE*qnorm(0.025)),2),
                         " – ",
                         round(exp(AoU_BETA + AoU_SE*qnorm(0.975)),2),
                         ")"),
         P = format.pval(P, digits = 2, eps = 1E-300),
         UKB_P = format.pval(P_UKB, digits = 2, eps = 1E-300),
         AoU_P = format.pval(P_AoU, digits = 2, eps = 1E-300),
         MAF_UKB_AoU = paste0(round(MAF_UKB,5)," | ",round(MAF_AoU,5))) %>%
  arrange(CHR, POS, Gene, Consequence, Cancer) %>%
  select(c(ID,Consequence,Gene,Cancer,Meta_OR,P,UKB_OR,UKB_P,AoU_OR,AoU_P,MAF_UKB_AoU))


write_table_for_word <- function(df, base_path = "table_output") {
  # ---- 1. Write standard TSV file ----
  tsv_file <- paste0(base_path, ".tsv")
  write.table(df, file = tsv_file, sep = "\t", quote = FALSE, row.names = FALSE)
  
  # ---- 2. Create Word-ready flextable ----
  ft <- flextable::flextable(df) %>%
    flextable::theme_box() %>%
    flextable::set_table_properties(layout = "autofit") %>%
    flextable::bold(part = "header") %>%
    flextable::fontsize(size = 6, part = "all") %>%
    flextable::autofit()
  
  # ---- 3. Save as Word document ----
  docx_file <- paste0(base_path, ".docx")
  flextable::save_as_docx(ft, path = docx_file)
  
}

write_table_for_word(formatted_sig_group_single_results,"austin_working/Dissertation/Meta/Aim1/scripts/result_summaries/tables/top_single_variant_individual_cancer_associations")


