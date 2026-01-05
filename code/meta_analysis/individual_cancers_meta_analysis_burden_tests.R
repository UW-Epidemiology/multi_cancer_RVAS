# Setup #
library(pacman)

p_load(rlang, tidyverse, data.table,
       R.utils, metafor, meta)

# Run meta-analyses #

meta_genes = c("ATM", "CHEK2", "BRCA2", "BRCA1", "CDKN2A", "TP53", "HAL", "HOXB13", "MSH6",  
               "OCA2", "DNMT3A", "MLH1", "PALB2", "CPVL", "FLG", "MITF", "PPM1D", "PRDM7", "RTEL1",  
               "SRSF2", "ASXL1", "FAM111A", "IFIH1", "NF1", "POT1", "TET2", "PLEKHA4", "HLA-DPB1",  
               "MICA", "JAK2", "IDH2", "IGLL5","SAMHD1")


meta_results = NULL

for (i in c("anus", "bladder", "bone", "breast",
            "brain", "cervix", "colorectal",
            "endometrium", "esophagus", "HL",
            "kidney", "leukemias", "liver",
            "lung", "melanoma", "myeloma",
            "NHL", "neck", "oral",
            "ovary", "pancreas", "prostate",
            "stomach", "testis", "thyroid",
            "eye","non_melanoma")) {
  
  pLoF_meta = NULL
  missense_meta = NULL
  synonymous_meta = NULL
  pLoF_missense_meta = NULL
  missense_synonymous_meta = NULL
  all_meta = NULL
  
  UKBB_Ns = fread(paste0("/lindstroem/austin_working/Dissertation/UKBB/gene_based_results/individual_cancers_combined/",i,"_SAIGE_results.txt.singleAssoc.txt"), nrows = 1) %>%
    select(N_case,N_ctrl)
  AoU_Ns = fread(paste0("/lindstroem/austin_working/Dissertation/AoU/gene_based_results/individual_cancers_combined/",i,"_SAIGE_results.txt.singleAssoc.txt"), nrows = 1) %>%
    select(N_case,N_ctrl)
  Ns = data.frame(N_case = UKBB_Ns$N_case + AoU_Ns$N_case,
                  N_ctrl = UKBB_Ns$N_ctrl + AoU_Ns$N_ctrl) %>% 
                    mutate(N_total = N_case + N_ctrl )
  
  
  UKBB_results = fread(paste0("/lindstroem/austin_working/Dissertation/UKBB/gene_based_results/individual_cancers_combined/",i,"_SAIGE_results.txt")) %>%
    filter(Region %in% meta_genes, MAC_case > 0, MAC_control > 0)
  AoU_results = fread(paste0("/lindstroem/austin_working/Dissertation/AoU/gene_based_results/individual_cancers_combined/",i,"_SAIGE_results.txt")) %>%
    filter(Region %in% meta_genes, MAC_case > 0, MAC_control > 0)
  
  r1_pLoF = UKBB_results %>% 
    filter(Group == "pLoF" & max_MAF == 0.5)
  r2_pLoF = AoU_results %>% 
    filter(Group == "pLoF" & max_MAF == 0.5)
  genes_r1_pLoF = unique(r1_pLoF$Region)
  genes_r2_pLoF = unique(r2_pLoF$Region)
  common_genes_pLoF = intersect(genes_r1_pLoF,genes_r2_pLoF)
  
  r1_missense = UKBB_results %>% 
    filter(Group == "missense" & max_MAF == 0.5)
  r2_missense = AoU_results %>% 
    filter(Group == "missense" & max_MAF == 0.5)
  genes_r1_missense = unique(r1_missense$Region)
  genes_r2_missense = unique(r2_missense$Region)
  common_genes_missense = intersect(genes_r1_missense,genes_r2_missense)
  
  r1_synonymous = UKBB_results %>% 
    filter(Group == "synonymous" & max_MAF == 0.01)
  r2_synonymous = AoU_results %>% 
    filter(Group == "synonymous" & max_MAF == 0.01)
  genes_r1_synonymous = unique(r1_synonymous$Region)
  genes_r2_synonymous = unique(r2_synonymous$Region)
  common_genes_synonymous = intersect(genes_r1_synonymous,genes_r2_synonymous)
  
  r1_pLoF_missense = UKBB_results %>% 
    filter(Group == "pLoF;missense" & max_MAF == 0.5)
  r2_pLoF_missense = AoU_results %>% 
    filter(Group == "pLoF;missense" & max_MAF == 0.5)
  genes_r1_pLoF_missense = unique(r1_pLoF_missense$Region)
  genes_r2_pLoF_missense = unique(r2_pLoF_missense$Region)
  common_genes_pLoF_missense = intersect(genes_r1_pLoF_missense,genes_r2_pLoF_missense)
  
  r1_missense_synonymous = UKBB_results %>% 
    filter(Group == "missense;synonymous" & max_MAF == 0.01)
  r2_missense_synonymous = AoU_results %>% 
    filter(Group == "missense;synonymous" & max_MAF == 0.01)
  genes_r1_missense_synonymous = unique(r1_missense_synonymous$Region)
  genes_r2_missense_synonymous = unique(r2_missense_synonymous$Region)
  common_genes_missense_synonymous = intersect(genes_r1_missense_synonymous,genes_r2_missense_synonymous)
  
  r1_all = UKBB_results %>% 
    filter(Group == "pLoF;missense;synonymous;other" & max_MAF == 0.01)
  r2_all = AoU_results %>% 
    filter(Group == "pLoF;missense;synonymous;other" & max_MAF == 0.01)
  genes_r1_all = unique(r1_all$Region)
  genes_r2_all = unique(r2_all$Region)
  common_genes_all = intersect(genes_r1_all,genes_r2_all)
  
  
  cancer_meta_pLoF_results = data.frame(Gene = NA,
                                Group = NA,
                                Cancer = NA,
                                N_case_meta = rep(Ns$N_case,length(common_genes_pLoF)),
                                N_ctrl_meta = rep(Ns$N_ctrl,length(common_genes_pLoF)),
                                N_total_meta = rep(Ns$N_total,length(common_genes_pLoF)),
                                BETA = NA,
                                BETA_SE = NA,
                                BETA_CI_LB = NA,
                                BETA_CI_UP = NA,
                                P = NA,
                                I2 = NA,
                                UKB_BETA = NA,
                                UKB_SE = NA,
                                UKB_P = NA,
                                AoU_BETA = NA,
                                AoU_SE = NA,
                                AoU_P = NA)

  for (j in 1:length(common_genes_pLoF)) {
    r1_pLoF_gene = r1_pLoF %>%
      filter(Region == common_genes_pLoF[j])
    r2_pLoF_gene = r2_pLoF %>%
      filter(Region == common_genes_pLoF[j])
    
    
    gene_meta_fe = rma(yi = c(r1_pLoF_gene$BETA_Burden,r2_pLoF_gene$BETA_Burden),
                    sei = c(r1_pLoF_gene$SE_Burden,r2_pLoF_gene$SE_Burden),method = "FE")
    
    cancer_meta_pLoF_results[j,1] = common_genes_pLoF[j]
    cancer_meta_pLoF_results[j,2] = r1_pLoF_gene$Group
    cancer_meta_pLoF_results[j,3] = i
    cancer_meta_pLoF_results[j,7] = gene_meta_fe$beta
    cancer_meta_pLoF_results[j,8] = gene_meta_fe$se
    cancer_meta_pLoF_results[j,9] = gene_meta_fe$ci.lb
    cancer_meta_pLoF_results[j,10] = gene_meta_fe$ci.ub
    cancer_meta_pLoF_results[j,11] = gene_meta_fe$pval
    cancer_meta_pLoF_results[j,12] = gene_meta_fe$I2
    cancer_meta_pLoF_results[j,13] = r1_pLoF_gene$BETA_Burden
    cancer_meta_pLoF_results[j,14] = r1_pLoF_gene$SE_Burden
    cancer_meta_pLoF_results[j,15] = r1_pLoF_gene$Pvalue_Burden
    cancer_meta_pLoF_results[j,16] = r2_pLoF_gene$BETA_Burden
    cancer_meta_pLoF_results[j,17] = r2_pLoF_gene$SE_Burden
    cancer_meta_pLoF_results[j,18] = r2_pLoF_gene$Pvalue_Burden
  }
  pLoF_meta = rbind(pLoF_meta,cancer_meta_pLoF_results)
  
  cancer_meta_missense_results = data.frame(Gene = NA,
                                     Group = NA,
                                     Cancer = NA,
                                     N_case_meta = rep(Ns$N_case,length(common_genes_missense)),
                                     N_ctrl_meta = rep(Ns$N_ctrl,length(common_genes_missense)),
                                     N_total_meta = rep(Ns$N_total,length(common_genes_missense)),
                                     BETA = NA,
                                     BETA_SE = NA,
                                     BETA_CI_LB = NA,
                                     BETA_CI_UP = NA,
                                     P = NA,
                                     I2 = NA,
                                     UKB_BETA = NA,
                                     UKB_SE = NA,
                                     UKB_P = NA,
                                     AoU_BETA = NA,
                                     AoU_SE = NA,
                                     AoU_P = NA)
  
  for (j in 1:length(common_genes_missense)) {
    r1_missense_gene = r1_missense %>%
      filter(Region == common_genes_missense[j])
    r2_missense_gene = r2_missense %>%
      filter(Region == common_genes_missense[j])
    
    
    gene_meta_fe = rma(yi = c(r1_missense_gene$BETA_Burden,r2_missense_gene$BETA_Burden),
                       sei = c(r1_missense_gene$SE_Burden,r2_missense_gene$SE_Burden),method = "FE")
    
    cancer_meta_missense_results[j,1] = common_genes_missense[j]
    cancer_meta_missense_results[j,2] = r1_missense_gene$Group
    cancer_meta_missense_results[j,3] = i
    cancer_meta_missense_results[j,7] = gene_meta_fe$beta
    cancer_meta_missense_results[j,8] = gene_meta_fe$se
    cancer_meta_missense_results[j,9] = gene_meta_fe$ci.lb
    cancer_meta_missense_results[j,10] = gene_meta_fe$ci.ub
    cancer_meta_missense_results[j,11] = gene_meta_fe$pval
    cancer_meta_missense_results[j,12] = gene_meta_fe$I2
    cancer_meta_missense_results[j,13] = r1_missense_gene$BETA_Burden
    cancer_meta_missense_results[j,14] = r1_missense_gene$SE_Burden
    cancer_meta_missense_results[j,15] = r1_missense_gene$Pvalue_Burden
    cancer_meta_missense_results[j,16] = r2_missense_gene$BETA_Burden
    cancer_meta_missense_results[j,17] = r2_missense_gene$SE_Burden
    cancer_meta_missense_results[j,18] = r2_missense_gene$Pvalue_Burden
  }
  missense_meta = rbind(missense_meta,cancer_meta_missense_results)
  
  cancer_meta_pLoF_missense_results = data.frame(Gene = NA,
                                         Group = NA,
                                         Cancer = NA,
                                         N_case_meta = rep(Ns$N_case,length(common_genes_pLoF_missense)),
                                         N_ctrl_meta = rep(Ns$N_ctrl,length(common_genes_pLoF_missense)),
                                         N_total_meta = rep(Ns$N_total,length(common_genes_pLoF_missense)),
                                         BETA = NA,
                                         BETA_SE = NA,
                                         BETA_CI_LB = NA,
                                         BETA_CI_UP = NA,
                                         P = NA,
                                         I2 = NA,
                                         UKB_BETA = NA,
                                         UKB_SE = NA,
                                         UKB_P = NA,
                                         AoU_BETA = NA,
                                         AoU_SE = NA,
                                         AoU_P = NA)
  
  for (j in 1:length(common_genes_pLoF_missense)) {
    r1_pLoF_missense_gene = r1_pLoF_missense %>%
      filter(Region == common_genes_pLoF_missense[j])
    r2_pLoF_missense_gene = r2_pLoF_missense %>%
      filter(Region == common_genes_pLoF_missense[j])
    
    
    gene_meta_fe = rma(yi = c(r1_pLoF_missense_gene$BETA_Burden,r2_pLoF_missense_gene$BETA_Burden),
                       sei = c(r1_pLoF_missense_gene$SE_Burden,r2_pLoF_missense_gene$SE_Burden),method = "FE")
    
    cancer_meta_pLoF_missense_results[j,1] = common_genes_pLoF_missense[j]
    cancer_meta_pLoF_missense_results[j,2] = r1_pLoF_missense_gene$Group
    cancer_meta_pLoF_missense_results[j,3] = i
    cancer_meta_pLoF_missense_results[j,7] = gene_meta_fe$beta
    cancer_meta_pLoF_missense_results[j,8] = gene_meta_fe$se
    cancer_meta_pLoF_missense_results[j,9] = gene_meta_fe$ci.lb
    cancer_meta_pLoF_missense_results[j,10] = gene_meta_fe$ci.ub
    cancer_meta_pLoF_missense_results[j,11] = gene_meta_fe$pval
    cancer_meta_pLoF_missense_results[j,12] = gene_meta_fe$I2
    cancer_meta_pLoF_missense_results[j,13] = r1_pLoF_missense_gene$BETA_Burden
    cancer_meta_pLoF_missense_results[j,14] = r1_pLoF_missense_gene$SE_Burden
    cancer_meta_pLoF_missense_results[j,15] = r1_pLoF_missense_gene$Pvalue_Burden
    cancer_meta_pLoF_missense_results[j,16] = r2_pLoF_missense_gene$BETA_Burden
    cancer_meta_pLoF_missense_results[j,17] = r2_pLoF_missense_gene$SE_Burden
    cancer_meta_pLoF_missense_results[j,18] = r2_pLoF_missense_gene$Pvalue_Burden
  }
  pLoF_missense_meta = rbind(pLoF_missense_meta,cancer_meta_pLoF_missense_results)
  
  cancer_meta_missense_synonymous_results = data.frame(Gene = NA,
                                         Group = NA,
                                         Cancer = NA,
                                         N_case_meta = rep(Ns$N_case,length(common_genes_missense_synonymous)),
                                         N_ctrl_meta = rep(Ns$N_ctrl,length(common_genes_missense_synonymous)),
                                         N_total_meta = rep(Ns$N_total,length(common_genes_missense_synonymous)),
                                         BETA = NA,
                                         BETA_SE = NA,
                                         BETA_CI_LB = NA,
                                         BETA_CI_UP = NA,
                                         P = NA,
                                         I2 = NA,
                                         UKB_BETA = NA,
                                         UKB_SE = NA,
                                         UKB_P = NA,
                                         AoU_BETA = NA,
                                         AoU_SE = NA,
                                         AoU_P = NA)
  
  for (j in 1:length(common_genes_missense_synonymous)) {
    r1_missense_synonymous_gene = r1_missense_synonymous %>%
      filter(Region == common_genes_missense_synonymous[j])
    r2_missense_synonymous_gene = r2_missense_synonymous %>%
      filter(Region == common_genes_missense_synonymous[j])
    
    
    gene_meta_fe = rma(yi = c(r1_missense_synonymous_gene$BETA_Burden,r2_missense_synonymous_gene$BETA_Burden),
                       sei = c(r1_missense_synonymous_gene$SE_Burden,r2_missense_synonymous_gene$SE_Burden),method = "FE")
    
    cancer_meta_missense_synonymous_results[j,1] = common_genes_missense_synonymous[j]
    cancer_meta_missense_synonymous_results[j,2] = r1_missense_synonymous_gene$Group
    cancer_meta_missense_synonymous_results[j,3] = i
    cancer_meta_missense_synonymous_results[j,7] = gene_meta_fe$beta
    cancer_meta_missense_synonymous_results[j,8] = gene_meta_fe$se
    cancer_meta_missense_synonymous_results[j,9] = gene_meta_fe$ci.lb
    cancer_meta_missense_synonymous_results[j,10] = gene_meta_fe$ci.ub
    cancer_meta_missense_synonymous_results[j,11] = gene_meta_fe$pval
    cancer_meta_missense_synonymous_results[j,12] = gene_meta_fe$I2
    cancer_meta_missense_synonymous_results[j,13] = r1_missense_synonymous_gene$BETA_Burden
    cancer_meta_missense_synonymous_results[j,14] = r1_missense_synonymous_gene$SE_Burden
    cancer_meta_missense_synonymous_results[j,15] = r1_missense_synonymous_gene$Pvalue_Burden
    cancer_meta_missense_synonymous_results[j,16] = r2_missense_synonymous_gene$BETA_Burden
    cancer_meta_missense_synonymous_results[j,17] = r2_missense_synonymous_gene$SE_Burden
    cancer_meta_missense_synonymous_results[j,18] = r2_missense_synonymous_gene$Pvalue_Burden
  }
  missense_synonymous_meta = rbind(missense_synonymous_meta,cancer_meta_missense_synonymous_results)
  
  cancer_meta_synonymous_results = data.frame(Gene = NA,
                                         Group = NA,
                                         Cancer = NA,
                                         N_case_meta = rep(Ns$N_case,length(common_genes_synonymous)),
                                         N_ctrl_meta = rep(Ns$N_ctrl,length(common_genes_synonymous)),
                                         N_total_meta = rep(Ns$N_total,length(common_genes_synonymous)),
                                         BETA = NA,
                                         BETA_SE = NA,
                                         BETA_CI_LB = NA,
                                         BETA_CI_UP = NA,
                                         P = NA,
                                         I2 = NA,
                                         UKB_BETA = NA,
                                         UKB_SE = NA,
                                         UKB_P = NA,
                                         AoU_BETA = NA,
                                         AoU_SE = NA,
                                         AoU_P = NA)
  
  for (j in 1:length(common_genes_synonymous)) {
    r1_synonymous_gene = r1_synonymous %>%
      filter(Region == common_genes_synonymous[j])
    r2_synonymous_gene = r2_synonymous %>%
      filter(Region == common_genes_synonymous[j])
    
    
    gene_meta_fe = rma(yi = c(r1_synonymous_gene$BETA_Burden,r2_synonymous_gene$BETA_Burden),
                       sei = c(r1_synonymous_gene$SE_Burden,r2_synonymous_gene$SE_Burden),method = "FE")
    
    cancer_meta_synonymous_results[j,1] = common_genes_synonymous[j]
    cancer_meta_synonymous_results[j,2] = r1_synonymous_gene$Group
    cancer_meta_synonymous_results[j,3] = i
    cancer_meta_synonymous_results[j,7] = gene_meta_fe$beta
    cancer_meta_synonymous_results[j,8] = gene_meta_fe$se
    cancer_meta_synonymous_results[j,9] = gene_meta_fe$ci.lb
    cancer_meta_synonymous_results[j,10] = gene_meta_fe$ci.ub
    cancer_meta_synonymous_results[j,11] = gene_meta_fe$pval
    cancer_meta_synonymous_results[j,12] = gene_meta_fe$I2
    cancer_meta_synonymous_results[j,13] = r1_synonymous_gene$BETA_Burden
    cancer_meta_synonymous_results[j,14] = r1_synonymous_gene$SE_Burden
    cancer_meta_synonymous_results[j,15] = r1_synonymous_gene$Pvalue_Burden
    cancer_meta_synonymous_results[j,16] = r2_synonymous_gene$BETA_Burden
    cancer_meta_synonymous_results[j,17] = r2_synonymous_gene$SE_Burden
    cancer_meta_synonymous_results[j,18] = r2_synonymous_gene$Pvalue_Burden
  }
  synonymous_meta = rbind(synonymous_meta,cancer_meta_synonymous_results)
  
  cancer_meta_all_results = data.frame(Gene = NA,
                                           Group = NA,
                                           Cancer = NA,
                                           N_case_meta = rep(Ns$N_case,length(common_genes_all)),
                                           N_ctrl_meta = rep(Ns$N_ctrl,length(common_genes_all)),
                                           N_total_meta = rep(Ns$N_total,length(common_genes_all)),
                                           BETA = NA,
                                           BETA_SE = NA,
                                           BETA_CI_LB = NA,
                                           BETA_CI_UP = NA,
                                           P = NA,
                                           I2 = NA,
                                           UKB_BETA = NA,
                                           UKB_SE = NA,
                                           UKB_P = NA,
                                           AoU_BETA = NA,
                                           AoU_SE = NA,
                                           AoU_P = NA)
  
  for (j in 1:length(common_genes_all)) {
    r1_all_gene = r1_all %>%
      filter(Region == common_genes_all[j])
    r2_all_gene = r2_all %>%
      filter(Region == common_genes_all[j])
    
    
    gene_meta_fe = rma(yi = c(r1_all_gene$BETA_Burden,r2_all_gene$BETA_Burden),
                       sei = c(r1_all_gene$SE_Burden,r2_all_gene$SE_Burden),method = "FE")
    
    cancer_meta_all_results[j,1] = common_genes_all[j]
    cancer_meta_all_results[j,2] = r1_all_gene$Group
    cancer_meta_all_results[j,3] = i
    cancer_meta_all_results[j,7] = gene_meta_fe$beta
    cancer_meta_all_results[j,8] = gene_meta_fe$se
    cancer_meta_all_results[j,9] = gene_meta_fe$ci.lb
    cancer_meta_all_results[j,10] = gene_meta_fe$ci.ub
    cancer_meta_all_results[j,11] = gene_meta_fe$pval
    cancer_meta_all_results[j,12] = gene_meta_fe$I2
    cancer_meta_all_results[j,13] = r1_all_gene$BETA_Burden
    cancer_meta_all_results[j,14] = r1_all_gene$SE_Burden
    cancer_meta_all_results[j,15] = r1_all_gene$Pvalue_Burden
    cancer_meta_all_results[j,16] = r2_all_gene$BETA_Burden
    cancer_meta_all_results[j,17] = r2_all_gene$SE_Burden
    cancer_meta_all_results[j,18] = r2_all_gene$Pvalue_Burden
  }
  all_meta = rbind(all_meta,cancer_meta_all_results)
  
  meta_results = rbind(meta_results,
                       pLoF_meta,
                       missense_meta,
                       synonymous_meta,
                       pLoF_missense_meta,
                       missense_synonymous_meta,
                       all_meta)
}

root_path <- "/lindstroem/austin_working/Dissertation/Meta/Aim1/individual_cancers"  

fwrite(meta_results, "/lindstroem/austin_working/Dissertation/Meta/Aim1/individual_cancers/individual_cancers_meta.txt",
       sep = "\t",
       quote = FALSE, na = "NA")
