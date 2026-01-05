# Setup #

library(pacman)

p_load(tidyverse, data.table, openxlsx)

# Load #

sig_genes = c("ATM", "CHEK2", "BRCA2", "SAMHD1", "BRCA1", "CDKN2A", "TP53", "HAL", "HOXB13", "MSH6",  
              "OCA2", "DNMT3A", "MLH1", "PALB2", "CPVL", "FLG", "MITF", "PPM1D", "PRDM7", "RTEL1",  
              "SRSF2", "ASXL1", "FAM111A", "IFIH1", "NF1", "POT1", "TET2", "PLEKHA4", "HLA-DPB1",  
              "MICA", "JAK2", "IDH2", "IGLL5")

meta_sig_genes = c("ATM", "CHEK2", "BRCA2", "SAMHD1", "BRCA1", "CDKN2A", "TP53", "HAL", "HOXB13", "MSH6",  
              "OCA2", "DNMT3A", "MLH1", "PALB2", "CPVL", "FLG", "MITF", "PPM1D", "PRDM7", "RTEL1",  
              "SRSF2", "ASXL1", "FAM111A", "IFIH1", "NF1", "POT1", "TET2")

UKB_sig_genes = c("ATM", "CHEK2", "BRCA2", "SAMHD1", "BRCA1", "CDKN2A", "TP53", "HAL", "HOXB13", "MSH6",  
                  "OCA2", "MLH1", "PALB2", "CPVL", "FLG", "MITF", "PPM1D", "PRDM7", "RTEL1",  
                  "SRSF2", "FAM111A", "IFIH1", "NF1", "POT1", "TET2", "PLEKHA4", "HLA-DPB1",  
                  "MICA")

AoU_sig_genes = c("ATM", "CHEK2", "BRCA2", "SAMHD1", "BRCA1", "CDKN2A", "TP53", "MLH1", "PALB2",
                  "PPM1D", "PRDM7", "NF1", "TET2")

group_data = NULL
for (i in c("ectoderm","mesoderm","endoderm","hormone","infectious","smoking")) {
  for (j in c("pLoF","missense","synonymous","pLoF_missense","missense_synonymous","all")) {
    data = fread(paste0("/lindstroem/austin_working/Dissertation/Meta/Aim1/cancer_groups/",j,"_meta_",i,".txt")) %>%
      filter(Gene %in% sig_genes) %>%
      mutate(Cancer = i)
    group_data = rbind.data.frame(group_data,data)
  }
}

# Significant meta-analysis burden test result tables

## File output path
out_file <- "./austin_working/Dissertation/Meta/Aim1/scripts/result_summaries/tables/meta_sig_burden_tables.xlsx"

## Define cancer types and corresponding sheet names
cancer_types <- c("ectoderm", "mesoderm", "endoderm", "hormone", "smoking", "infectious")
sheet_names <- c("Ectoderm", "Mesoderm", "Endoderm", "Hormone", "Smoking", "Infectious")

## Create a new workbook
wb <- createWorkbook()

## Loop through each cancer type and add a worksheet
for (i in seq_along(cancer_types)) {
  cancer_type <- cancer_types[i]
  sheet_name <- sheet_names[i]
  
  table <- group_data %>%
    filter(Gene %in% meta_sig_genes, Cancer == cancer_type, P < 4.79E-7) %>%
    group_by(Gene, Group) %>%
    mutate(
      OR_CI = paste0(
        round(exp(BETA), 4), " (",
        round(exp(BETA_CI_LB), 4), ", ",
        round(exp(BETA_CI_UP), 4), ")"
      )
    ) %>%
    ungroup() %>%
    mutate(P = formatC(P, format = "E", digits = 2)) %>%
    arrange(Gene, Group) %>%
    select(Gene, Group, OR_CI, P)
  
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet = sheet_name, table)
}

## Save the workbook
saveWorkbook(wb, out_file, overwrite = TRUE)

# Significant UKB burden test result tables

## File output path
out_file <- "./austin_working/Dissertation/Meta/Aim1/scripts/result_summaries/tables/ukb_sig_burden_tables.xlsx"

## Define cancer types and corresponding sheet names
cancer_types <- c("ectoderm", "mesoderm", "endoderm", "hormone", "smoking", "infectious")
sheet_names <- c("Ectoderm", "Mesoderm", "Endoderm", "Hormone", "Smoking", "Infectious")

## Create a new workbook
wb <- createWorkbook()

## Loop through each cancer type and add a worksheet
for (i in seq_along(cancer_types)) {
  cancer_type <- cancer_types[i]
  sheet_name <- sheet_names[i]
  
  table <- group_data %>%
    filter(Gene %in% UKB_sig_genes,
           Cancer == cancer_type, UKB_P < 4.77E-7) %>%
    group_by(Gene, Group) %>%
    mutate(
      OR_CI = paste0(
        round(exp(UKB_BETA), 4), " (",
        round(exp(UKB_BETA+(UKB_SE*qnorm(0.025))), 4), ", ",
        round(exp(UKB_BETA+(UKB_SE*qnorm(0.975))), 4), ")"
      )
    ) %>%
    ungroup() %>%
    mutate(P = formatC(P, format = "E", digits = 2)) %>%
    arrange(Gene, Group) %>%
    select(Gene, Group, OR_CI, P)
  
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet = sheet_name, table)
}

## Save the workbook
saveWorkbook(wb, out_file, overwrite = TRUE)

# Significant AoU burden test result tables

## File output path
out_file <- "./austin_working/Dissertation/Meta/Aim1/scripts/result_summaries/tables/aou_sig_burden_tables.xlsx"

## Define cancer types and corresponding sheet names
cancer_types <- c("ectoderm", "mesoderm", "endoderm", "hormone", "smoking", "infectious")
sheet_names <- c("Ectoderm", "Mesoderm", "Endoderm", "Hormone", "Smoking", "Infectious")

## Create a new workbook
wb <- createWorkbook()

## Loop through each cancer type and add a worksheet
for (i in seq_along(cancer_types)) {
  cancer_type <- cancer_types[i]
  sheet_name <- sheet_names[i]
  
  table <- group_data %>%
    filter(Gene %in% AoU_sig_genes,
           Cancer == cancer_type, AoU_P < 4.68E-7) %>%
    group_by(Gene, Group) %>%
    mutate(
      OR_CI = paste0(
        round(exp(AoU_BETA), 4), " (",
        round(exp(AoU_BETA+(AoU_SE*qnorm(0.025))), 4), ", ",
        round(exp(AoU_BETA+(AoU_SE*qnorm(0.975))), 4), ")"
      )
    ) %>%
    ungroup() %>%
    mutate(P = formatC(P, format = "E", digits = 2)) %>%
    arrange(Gene, Group) %>%
    select(Gene, Group, OR_CI, P)
  
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet = sheet_name, table)
}

## Save the workbook
saveWorkbook(wb, out_file, overwrite = TRUE)


