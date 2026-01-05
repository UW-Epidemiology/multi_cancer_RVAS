# Setup #
library(pacman)

p_load(rlang, tidyverse, data.table,
       R.utils, metafor, meta)

# Run meta-analyses #

meta_genes = c("ATM", "CHEK2", "BRCA2", "BRCA1", "CDKN2A", "TP53", "HAL", "HOXB13", "MSH6",  
               "OCA2", "DNMT3A", "MLH1", "PALB2", "CPVL", "FLG", "MITF", "PPM1D", "PRDM7", "RTEL1",  
               "SRSF2", "ASXL1", "FAM111A", "IFIH1", "NF1", "POT1", "TET2", "PLEKHA4", "HLA-DPB1", "HLA.DPB1",
               "MICA", "JAK2", "IDH2", "IGLL5","SAMHD1")


## Run for regular carrier status analyses ##
meta_results = NULL
pLoF_meta = NULL
missense_meta = NULL


for (i in c("anus", "bladder", "bone", "breast",
            "brain", "cervix", "colorectal",
            "endometrium", "esophagus", "HL",
            "kidney", "leukemias", "liver",
            "lung", "melanoma", "myeloma",
            "NHL", "neck", "oral",
            "ovary", "pancreas", "prostate",
            "stomach", "testis", "thyroid",
            "eye","non_melanoma")) {
  
  UKBB_results =  NULL
  for (chr in c(1,2,3,4,6,7,9,11,12,13,15,16,17,19,20,22)) {
    chr_results = fread(paste0("/lindstroem/austin_working/Dissertation/UKBB/gene_based_results/individual_cancers_carrier_status/chr",chr,"_carrier_status_results.txt")) %>%
      filter(Gene %in% meta_genes, Cancer == i)
    UKBB_results = rbind.data.frame(UKBB_results,chr_results)
  }
  
  AoU_results = NULL
  for (chr in c(1,2,3,4,6,7,9,11,12,13,15,16,17,19,20,22)) {
    chr_results = fread(paste0("/lindstroem/austin_working/Dissertation/AoU/gene_based_results/individual_cancers_carrier_status/chr",chr,"_carrier_status_results.txt")) %>%
      filter(Gene %in% meta_genes, Cancer == i)
    chr_results$Gene <- gsub("\\.", "-", chr_results$Gene)
    AoU_results = rbind.data.frame(AoU_results, chr_results)
  }
  
  Ns = data.frame(N_case = median(UKBB_results$N_cases)+ median(AoU_results$N_cases),
                  N_ctrl = median(UKBB_results$N_controls) + median(AoU_results$N_controls)) %>% 
    mutate(N_total = N_case + N_ctrl)
  
  
  r1_pLoF = UKBB_results %>% 
    filter(Group == "pLoF")
  r2_pLoF = AoU_results %>% 
    filter(Group == "pLoF")
  genes_r1_pLoF = unique(r1_pLoF$Gene)
  genes_r2_pLoF = unique(r2_pLoF$Gene)
  common_genes_pLoF = intersect(genes_r1_pLoF,genes_r2_pLoF)
  
  r1_missense = UKBB_results %>% 
    filter(Group == "missense")
  r2_missense = AoU_results %>% 
    filter(Group == "missense")
  genes_r1_missense = unique(r1_missense$Gene)
  genes_r2_missense = unique(r2_missense$Gene)
  common_genes_missense = intersect(genes_r1_missense,genes_r2_missense)
  
  cancer_meta_pLoF_results = data.frame(Gene = NA,
                                Group = NA,
                                Cancer = NA,
                                N_case_meta = rep(Ns$N_case,length(common_genes_pLoF)),
                                N_ctrl_meta = rep(Ns$N_ctrl,length(common_genes_pLoF)),
                                N_total_meta = rep(Ns$N_total,length(common_genes_pLoF)),
                                Carrier_log_OR = NA,
                                Carrier_log_OR_SE = NA,
                                Carrier_OR = NA,
                                Carrier_OR_CI_LB = NA,
                                Carrier_OR_CI_UP = NA,
                                Carrier_P = NA,
                                I2 = NA,
                                Carrier_UKB_OR = NA,
                                Carrier_UKB_CI_LB = NA,
                                Carrier_UKB_CI_UP = NA,
                                Carrier_UKB_P = NA,
                                Carrier_AoU_OR = NA,
                                Carrier_AoU_CI_LB = NA,
                                Carrier_AoU_CI_UP = NA,
                                Carrier_AoU_P = NA)

  for (j in 1:length(common_genes_pLoF)) {
    r1_pLoF_gene = r1_pLoF %>%
      filter(Gene == common_genes_pLoF[j])
    r2_pLoF_gene = r2_pLoF %>%
      filter(Gene == common_genes_pLoF[j])
    
    
    gene_meta_fe = rma(yi = c(r1_pLoF_gene$log_OR,r2_pLoF_gene$log_OR),
                    sei = c(r1_pLoF_gene$log_OR_SE,r2_pLoF_gene$log_OR_SE),method = "FE")
    
    cancer_meta_pLoF_results[j,1] = common_genes_pLoF[j]
    cancer_meta_pLoF_results[j,2] = r1_pLoF_gene$Group
    cancer_meta_pLoF_results[j,3] = i
    cancer_meta_pLoF_results[j,7] = gene_meta_fe$beta
    cancer_meta_pLoF_results[j,8] = gene_meta_fe$se
    cancer_meta_pLoF_results[j,9] = exp(gene_meta_fe$beta)
    cancer_meta_pLoF_results[j,10] = exp(gene_meta_fe$ci.lb)
    cancer_meta_pLoF_results[j,11] = exp(gene_meta_fe$ci.ub)
    cancer_meta_pLoF_results[j,12] = gene_meta_fe$pval
    cancer_meta_pLoF_results[j,13] = gene_meta_fe$I2
    cancer_meta_pLoF_results[j,14] = r1_pLoF_gene$OR
    cancer_meta_pLoF_results[j,15] = r1_pLoF_gene$CI_lower
    cancer_meta_pLoF_results[j,16] = r1_pLoF_gene$CI_upper
    cancer_meta_pLoF_results[j,17] = r1_pLoF_gene$P
    cancer_meta_pLoF_results[j,18] = r2_pLoF_gene$OR
    cancer_meta_pLoF_results[j,19] = r2_pLoF_gene$CI_lower
    cancer_meta_pLoF_results[j,20] = r2_pLoF_gene$CI_upper
    cancer_meta_pLoF_results[j,21] = r2_pLoF_gene$P
  }
  pLoF_meta = rbind(pLoF_meta,cancer_meta_pLoF_results)
  
  cancer_meta_missense_results = data.frame(Gene = NA,
                                        Group = NA,
                                        Cancer = NA,
                                        N_case_meta = rep(Ns$N_case,length(common_genes_missense)),
                                        N_ctrl_meta = rep(Ns$N_ctrl,length(common_genes_missense)),
                                        N_total_meta = rep(Ns$N_total,length(common_genes_missense)),
                                        Carrier_log_OR = NA,
                                        Carrier_log_OR_SE = NA,
                                        Carrier_OR = NA,
                                        Carrier_OR_CI_LB = NA,
                                        Carrier_OR_CI_UP = NA,
                                        Carrier_P = NA,
                                        I2 = NA,
                                        Carrier_UKB_OR = NA,
                                        Carrier_UKB_CI_LB = NA,
                                        Carrier_UKB_CI_UP = NA,
                                        Carrier_UKB_P = NA,
                                        Carrier_AoU_OR = NA,
                                        Carrier_AoU_CI_LB = NA,
                                        Carrier_AoU_CI_UP = NA,
                                        Carrier_AoU_P = NA)  
  for (j in 1:length(common_genes_missense)) {
    r1_missense_gene = r1_missense %>%
      filter(Gene == common_genes_missense[j])
    r2_missense_gene = r2_missense %>%
      filter(Gene == common_genes_missense[j])
    
    
    gene_meta_fe = rma(yi = c(r1_missense_gene$log_OR,r2_missense_gene$log_OR),
                       sei = c(r1_missense_gene$log_OR_SE,r2_missense_gene$log_OR_SE),method = "FE")
    
    cancer_meta_missense_results[j,1] = common_genes_missense[j]
    cancer_meta_missense_results[j,2] = r1_missense_gene$Group
    cancer_meta_missense_results[j,3] = i
    cancer_meta_missense_results[j,7] = gene_meta_fe$beta
    cancer_meta_missense_results[j,8] = gene_meta_fe$se
    cancer_meta_missense_results[j,9] = exp(gene_meta_fe$beta)
    cancer_meta_missense_results[j,10] = exp(gene_meta_fe$ci.lb)
    cancer_meta_missense_results[j,11] = exp(gene_meta_fe$ci.ub)
    cancer_meta_missense_results[j,12] = gene_meta_fe$pval
    cancer_meta_missense_results[j,13] = gene_meta_fe$I2
    cancer_meta_missense_results[j,14] = r1_missense_gene$OR
    cancer_meta_missense_results[j,15] = r1_missense_gene$CI_lower
    cancer_meta_missense_results[j,16] = r1_missense_gene$CI_upper
    cancer_meta_missense_results[j,17] = r1_missense_gene$P
    cancer_meta_missense_results[j,18] = r2_missense_gene$OR
    cancer_meta_missense_results[j,19] = r2_missense_gene$CI_lower
    cancer_meta_missense_results[j,20] = r2_missense_gene$CI_upper
    cancer_meta_missense_results[j,21] = r2_missense_gene$P
  }
  missense_meta = rbind(missense_meta,cancer_meta_missense_results)
  meta_results = rbind.data.frame(pLoF_meta, missense_meta)
}

  
fwrite(meta_results, "/lindstroem/austin_working/Dissertation/Meta/Aim1/individual_cancers/individual_cancers_carrier_status_meta.txt",
       sep = "\t",
       quote = FALSE, na = "NA")


## Run for regular ClinVar carrier status analyses ##

meta_results = NULL
pLoF_meta = NULL
missense_meta = NULL


for (i in c("anus", "bladder", "bone", "breast",
            "brain", "cervix", "colorectal",
            "endometrium", "esophagus", "HL",
            "kidney", "leukemias", "liver",
            "lung", "melanoma", "myeloma",
            "NHL", "neck", "oral",
            "ovary", "pancreas", "prostate",
            "stomach", "testis", "thyroid",
            "eye","non_melanoma")) {
  
  UKBB_results =  NULL
  for (chr in c(1,2,3,4,7,9,11,13,15,16,17,20,22)) {
    chr_results = fread(paste0("/lindstroem/austin_working/Dissertation/UKBB/gene_based_results/individual_cancers_carrier_status/ClinVar_results/chr",chr,"_ClinVar_carrier_status_results.txt")) %>%
      filter(Gene %in% meta_genes, Cancer == i)
    UKBB_results = rbind.data.frame(UKBB_results,chr_results)
  }
  
  AoU_results = NULL
  for (chr in c(1,2,3,4,7,9,11,13,15,16,17,20,22)) {
    chr_results = fread(paste0("/lindstroem/austin_working/Dissertation/AoU/gene_based_results/individual_cancers_carrier_status//ClinVar_results/chr",chr,"_ClinVar_carrier_status_results.txt")) %>%
      filter(Gene %in% meta_genes, Cancer == i)
    chr_results$Gene <- gsub("\\.", "-", chr_results$Gene)
    AoU_results = rbind.data.frame(AoU_results,chr_results)
  }
  
  Ns = data.frame(N_case = median(UKBB_results$N_cases)+ median(AoU_results$N_cases),
                  N_ctrl = median(UKBB_results$N_controls) + median(AoU_results$N_controls)) %>% 
    mutate(N_total = N_case + N_ctrl)
  
  
  r1_pLoF = UKBB_results %>% 
    filter(Group == "pLoF")
  r2_pLoF = AoU_results %>% 
    filter(Group == "pLoF")
  genes_r1_pLoF = unique(r1_pLoF$Gene)
  genes_r2_pLoF = unique(r2_pLoF$Gene)
  common_genes_pLoF = intersect(genes_r1_pLoF,genes_r2_pLoF)
  
  r1_missense = UKBB_results %>% 
    filter(Group == "missense")
  r2_missense = AoU_results %>% 
    filter(Group == "missense")
  genes_r1_missense = unique(r1_missense$Gene)
  genes_r2_missense = unique(r2_missense$Gene)
  common_genes_missense = intersect(genes_r1_missense,genes_r2_missense)
  
  cancer_meta_pLoF_results = data.frame(Gene = NA,
                                        Group = NA,
                                        Cancer = NA,
                                        N_case_meta = rep(Ns$N_case,length(common_genes_pLoF)),
                                        N_ctrl_meta = rep(Ns$N_ctrl,length(common_genes_pLoF)),
                                        N_total_meta = rep(Ns$N_total,length(common_genes_pLoF)),
                                        ClinVar_Carrier_log_OR = NA,
                                        ClinVar_Carrier_log_OR_SE = NA,
                                        ClinVar_Carrier_OR = NA,
                                        ClinVar_Carrier_OR_CI_LB = NA,
                                        ClinVar_Carrier_OR_CI_UP = NA,
                                        ClinVar_Carrier_P = NA,
                                        I2 = NA,
                                        ClinVar_Carrier_UKB_OR = NA,
                                        ClinVar_Carrier_UKB_CI_LB = NA,
                                        ClinVar_Carrier_UKB_CI_UP = NA,
                                        ClinVar_Carrier_UKB_P = NA,
                                        ClinVar_Carrier_AoU_OR = NA,
                                        ClinVar_Carrier_AoU_CI_LB = NA,
                                        ClinVar_Carrier_AoU_CI_UP = NA,
                                        ClinVar_Carrier_AoU_P = NA)
  
  for (j in 1:length(common_genes_pLoF)) {
    r1_pLoF_gene = r1_pLoF %>%
      filter(Gene == common_genes_pLoF[j])
    r2_pLoF_gene = r2_pLoF %>%
      filter(Gene == common_genes_pLoF[j])
    
    
    gene_meta_fe = rma(yi = c(r1_pLoF_gene$log_OR,r2_pLoF_gene$log_OR),
                       sei = c(r1_pLoF_gene$log_OR_SE,r2_pLoF_gene$log_OR_SE),method = "FE")
    
    cancer_meta_pLoF_results[j,1] = common_genes_pLoF[j]
    cancer_meta_pLoF_results[j,2] = r1_pLoF_gene$Group
    cancer_meta_pLoF_results[j,3] = i
    cancer_meta_pLoF_results[j,7] = gene_meta_fe$beta
    cancer_meta_pLoF_results[j,8] = gene_meta_fe$se
    cancer_meta_pLoF_results[j,9] = exp(gene_meta_fe$beta)
    cancer_meta_pLoF_results[j,10] = exp(gene_meta_fe$ci.lb)
    cancer_meta_pLoF_results[j,11] = exp(gene_meta_fe$ci.ub)
    cancer_meta_pLoF_results[j,12] = gene_meta_fe$pval
    cancer_meta_pLoF_results[j,13] = gene_meta_fe$I2
    cancer_meta_pLoF_results[j,14] = r1_pLoF_gene$OR
    cancer_meta_pLoF_results[j,15] = r1_pLoF_gene$CI_lower
    cancer_meta_pLoF_results[j,16] = r1_pLoF_gene$CI_upper
    cancer_meta_pLoF_results[j,17] = r1_pLoF_gene$P
    cancer_meta_pLoF_results[j,18] = r2_pLoF_gene$OR
    cancer_meta_pLoF_results[j,19] = r2_pLoF_gene$CI_lower
    cancer_meta_pLoF_results[j,20] = r2_pLoF_gene$CI_upper
    cancer_meta_pLoF_results[j,21] = r2_pLoF_gene$P
  }
  pLoF_meta = rbind(pLoF_meta,cancer_meta_pLoF_results)
  
  cancer_meta_missense_results = data.frame(Gene = NA,
                                            Group = NA,
                                            Cancer = NA,
                                            N_case_meta = rep(Ns$N_case,length(common_genes_missense)),
                                            N_ctrl_meta = rep(Ns$N_ctrl,length(common_genes_missense)),
                                            N_total_meta = rep(Ns$N_total,length(common_genes_missense)),
                                            ClinVar_Carrier_log_OR = NA,
                                            ClinVar_Carrier_log_OR_SE = NA,
                                            ClinVar_Carrier_OR = NA,
                                            ClinVar_Carrier_OR_CI_LB = NA,
                                            ClinVar_Carrier_OR_CI_UP = NA,
                                            ClinVar_Carrier_P = NA,
                                            I2 = NA,
                                            ClinVar_Carrier_UKB_OR = NA,
                                            ClinVar_Carrier_UKB_CI_LB = NA,
                                            ClinVar_Carrier_UKB_CI_UP = NA,
                                            ClinVar_Carrier_UKB_P = NA,
                                            ClinVar_Carrier_AoU_OR = NA,
                                            ClinVar_Carrier_AoU_CI_LB = NA,
                                            ClinVar_Carrier_AoU_CI_UP = NA,
                                            ClinVar_Carrier_AoU_P = NA)  
  for (j in 1:length(common_genes_missense)) {
    r1_missense_gene = r1_missense %>%
      filter(Gene == common_genes_missense[j])
    r2_missense_gene = r2_missense %>%
      filter(Gene == common_genes_missense[j])
    
    
    gene_meta_fe = rma(yi = c(r1_missense_gene$log_OR,r2_missense_gene$log_OR),
                       sei = c(r1_missense_gene$log_OR_SE,r2_missense_gene$log_OR_SE),method = "FE")
    
    cancer_meta_missense_results[j,1] = common_genes_missense[j]
    cancer_meta_missense_results[j,2] = r1_missense_gene$Group
    cancer_meta_missense_results[j,3] = i
    cancer_meta_missense_results[j,7] = gene_meta_fe$beta
    cancer_meta_missense_results[j,8] = gene_meta_fe$se
    cancer_meta_missense_results[j,9] = exp(gene_meta_fe$beta)
    cancer_meta_missense_results[j,10] = exp(gene_meta_fe$ci.lb)
    cancer_meta_missense_results[j,11] = exp(gene_meta_fe$ci.ub)
    cancer_meta_missense_results[j,12] = gene_meta_fe$pval
    cancer_meta_missense_results[j,13] = gene_meta_fe$I2
    cancer_meta_missense_results[j,14] = r1_missense_gene$OR
    cancer_meta_missense_results[j,15] = r1_missense_gene$CI_lower
    cancer_meta_missense_results[j,16] = r1_missense_gene$CI_upper
    cancer_meta_missense_results[j,17] = r1_missense_gene$P
    cancer_meta_missense_results[j,18] = r2_missense_gene$OR
    cancer_meta_missense_results[j,19] = r2_missense_gene$CI_lower
    cancer_meta_missense_results[j,20] = r2_missense_gene$CI_upper
    cancer_meta_missense_results[j,21] = r2_missense_gene$P
  }
  missense_meta = rbind(missense_meta,cancer_meta_missense_results)
  meta_results = rbind.data.frame(pLoF_meta, missense_meta)
}


fwrite(meta_results, "/lindstroem/austin_working/Dissertation/Meta/Aim1/individual_cancers/individual_cancers_ClinVar_carrier_status_meta.txt",
       sep = "\t",
       quote = FALSE, na = "NA")
