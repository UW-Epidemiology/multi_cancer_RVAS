# Setup #
library(pacman)

p_load(rlang, tidyverse, data.table,
       R.utils, metafor, meta)

# Run ectoderm meta-analyses #

pLoF_meta = NULL
missense_meta = NULL
synonymous_meta = NULL
pLoF_missense_meta = NULL
missense_synonymous_meta = NULL
all_meta = NULL

for (i in 1:22) {
  UKBB_Ns = fread("/lindstroem/austin_working/Dissertation/UKBB/gene_based_results/cancer_groups/sparse_GRM/chr21_ectoderm_SAIGE_results.txt.singleAssoc.txt", nrows = 1) %>%
    select(N_case,N_ctrl)
  AoU_Ns = fread("/lindstroem/austin_working/Dissertation/AoU/gene_based_results/cancer_groups/sparse_GRM/chr21_ectoderm_SAIGE_results.txt.singleAssoc.txt", nrows = 1) %>%
    select(N_case,N_ctrl)
  Ns = data.frame(N_case = UKBB_Ns$N_case + AoU_Ns$N_case,
                  N_ctrl = UKBB_Ns$N_ctrl + AoU_Ns$N_ctrl) %>% 
                    mutate(N_total = N_case + N_ctrl )
  
  
  UKBB_results = fread(paste0("/lindstroem/austin_working/Dissertation/UKBB/gene_based_results/cancer_groups/sparse_GRM/chr",i,"_ectoderm_SAIGE_results.txt"))
  AoU_results = fread(paste0("/lindstroem/austin_working/Dissertation/AoU/gene_based_results/cancer_groups/sparse_GRM/chr",i,"_ectoderm_SAIGE_results.txt"))
  
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
  
  
  chr_meta_pLoF_results = data.frame(Gene = NA,
                                Group = NA,
                                CHR = NA,
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
    
    chr_meta_pLoF_results[j,1] = common_genes_pLoF[j]
    chr_meta_pLoF_results[j,2] = r1_pLoF_gene$Group
    chr_meta_pLoF_results[j,3] = i
    chr_meta_pLoF_results[j,7] = gene_meta_fe$beta
    chr_meta_pLoF_results[j,8] = gene_meta_fe$se
    chr_meta_pLoF_results[j,9] = gene_meta_fe$ci.lb
    chr_meta_pLoF_results[j,10] = gene_meta_fe$ci.ub
    chr_meta_pLoF_results[j,11] = gene_meta_fe$pval
    chr_meta_pLoF_results[j,12] = gene_meta_fe$I2
    chr_meta_pLoF_results[j,13] = r1_pLoF_gene$BETA_Burden
    chr_meta_pLoF_results[j,14] = r1_pLoF_gene$SE_Burden
    chr_meta_pLoF_results[j,15] = r1_pLoF_gene$Pvalue_Burden
    chr_meta_pLoF_results[j,16] = r2_pLoF_gene$BETA_Burden
    chr_meta_pLoF_results[j,17] = r2_pLoF_gene$SE_Burden
    chr_meta_pLoF_results[j,18] = r2_pLoF_gene$Pvalue_Burden
  }
  pLoF_meta = rbind(pLoF_meta,chr_meta_pLoF_results)
  
  chr_meta_missense_results = data.frame(Gene = NA,
                                     Group = NA,
                                     CHR = NA,
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
    
    chr_meta_missense_results[j,1] = common_genes_missense[j]
    chr_meta_missense_results[j,2] = r1_missense_gene$Group
    chr_meta_missense_results[j,3] = i
    chr_meta_missense_results[j,7] = gene_meta_fe$beta
    chr_meta_missense_results[j,8] = gene_meta_fe$se
    chr_meta_missense_results[j,9] = gene_meta_fe$ci.lb
    chr_meta_missense_results[j,10] = gene_meta_fe$ci.ub
    chr_meta_missense_results[j,11] = gene_meta_fe$pval
    chr_meta_missense_results[j,12] = gene_meta_fe$I2
    chr_meta_missense_results[j,13] = r1_missense_gene$BETA_Burden
    chr_meta_missense_results[j,14] = r1_missense_gene$SE_Burden
    chr_meta_missense_results[j,15] = r1_missense_gene$Pvalue_Burden
    chr_meta_missense_results[j,16] = r2_missense_gene$BETA_Burden
    chr_meta_missense_results[j,17] = r2_missense_gene$SE_Burden
    chr_meta_missense_results[j,18] = r2_missense_gene$Pvalue_Burden
  }
  missense_meta = rbind(missense_meta,chr_meta_missense_results)
  
  chr_meta_pLoF_missense_results = data.frame(Gene = NA,
                                         Group = NA,
                                         CHR = NA,
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
    
    chr_meta_pLoF_missense_results[j,1] = common_genes_pLoF_missense[j]
    chr_meta_pLoF_missense_results[j,2] = r1_pLoF_missense_gene$Group
    chr_meta_pLoF_missense_results[j,3] = i
    chr_meta_pLoF_missense_results[j,7] = gene_meta_fe$beta
    chr_meta_pLoF_missense_results[j,8] = gene_meta_fe$se
    chr_meta_pLoF_missense_results[j,9] = gene_meta_fe$ci.lb
    chr_meta_pLoF_missense_results[j,10] = gene_meta_fe$ci.ub
    chr_meta_pLoF_missense_results[j,11] = gene_meta_fe$pval
    chr_meta_pLoF_missense_results[j,12] = gene_meta_fe$I2
    chr_meta_pLoF_missense_results[j,13] = r1_pLoF_missense_gene$BETA_Burden
    chr_meta_pLoF_missense_results[j,14] = r1_pLoF_missense_gene$SE_Burden
    chr_meta_pLoF_missense_results[j,15] = r1_pLoF_missense_gene$Pvalue_Burden
    chr_meta_pLoF_missense_results[j,16] = r2_pLoF_missense_gene$BETA_Burden
    chr_meta_pLoF_missense_results[j,17] = r2_pLoF_missense_gene$SE_Burden
    chr_meta_pLoF_missense_results[j,18] = r2_pLoF_missense_gene$Pvalue_Burden
  }
  pLoF_missense_meta = rbind(pLoF_missense_meta,chr_meta_pLoF_missense_results)
  
  chr_meta_missense_synonymous_results = data.frame(Gene = NA,
                                         Group = NA,
                                         CHR = NA,
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
    
    chr_meta_missense_synonymous_results[j,1] = common_genes_missense_synonymous[j]
    chr_meta_missense_synonymous_results[j,2] = r1_missense_synonymous_gene$Group
    chr_meta_missense_synonymous_results[j,3] = i
    chr_meta_missense_synonymous_results[j,7] = gene_meta_fe$beta
    chr_meta_missense_synonymous_results[j,8] = gene_meta_fe$se
    chr_meta_missense_synonymous_results[j,9] = gene_meta_fe$ci.lb
    chr_meta_missense_synonymous_results[j,10] = gene_meta_fe$ci.ub
    chr_meta_missense_synonymous_results[j,11] = gene_meta_fe$pval
    chr_meta_missense_synonymous_results[j,12] = gene_meta_fe$I2
    chr_meta_missense_synonymous_results[j,13] = r1_missense_synonymous_gene$BETA_Burden
    chr_meta_missense_synonymous_results[j,14] = r1_missense_synonymous_gene$SE_Burden
    chr_meta_missense_synonymous_results[j,15] = r1_missense_synonymous_gene$Pvalue_Burden
    chr_meta_missense_synonymous_results[j,16] = r2_missense_synonymous_gene$BETA_Burden
    chr_meta_missense_synonymous_results[j,17] = r2_missense_synonymous_gene$SE_Burden
    chr_meta_missense_synonymous_results[j,18] = r2_missense_synonymous_gene$Pvalue_Burden
  }
  missense_synonymous_meta = rbind(missense_synonymous_meta,chr_meta_missense_synonymous_results)
  
  chr_meta_synonymous_results = data.frame(Gene = NA,
                                         Group = NA,
                                         CHR = NA,
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
    
    chr_meta_synonymous_results[j,1] = common_genes_synonymous[j]
    chr_meta_synonymous_results[j,2] = r1_synonymous_gene$Group
    chr_meta_synonymous_results[j,3] = i
    chr_meta_synonymous_results[j,7] = gene_meta_fe$beta
    chr_meta_synonymous_results[j,8] = gene_meta_fe$se
    chr_meta_synonymous_results[j,9] = gene_meta_fe$ci.lb
    chr_meta_synonymous_results[j,10] = gene_meta_fe$ci.ub
    chr_meta_synonymous_results[j,11] = gene_meta_fe$pval
    chr_meta_synonymous_results[j,12] = gene_meta_fe$I2
    chr_meta_synonymous_results[j,13] = r1_synonymous_gene$BETA_Burden
    chr_meta_synonymous_results[j,14] = r1_synonymous_gene$SE_Burden
    chr_meta_synonymous_results[j,15] = r1_synonymous_gene$Pvalue_Burden
    chr_meta_synonymous_results[j,16] = r2_synonymous_gene$BETA_Burden
    chr_meta_synonymous_results[j,17] = r2_synonymous_gene$SE_Burden
    chr_meta_synonymous_results[j,18] = r2_synonymous_gene$Pvalue_Burden
  }
  synonymous_meta = rbind(synonymous_meta,chr_meta_synonymous_results)
  
  chr_meta_all_results = data.frame(Gene = NA,
                                           Group = NA,
                                           CHR = NA,
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
    
    chr_meta_all_results[j,1] = common_genes_all[j]
    chr_meta_all_results[j,2] = r1_all_gene$Group
    chr_meta_all_results[j,3] = i
    chr_meta_all_results[j,7] = gene_meta_fe$beta
    chr_meta_all_results[j,8] = gene_meta_fe$se
    chr_meta_all_results[j,9] = gene_meta_fe$ci.lb
    chr_meta_all_results[j,10] = gene_meta_fe$ci.ub
    chr_meta_all_results[j,11] = gene_meta_fe$pval
    chr_meta_all_results[j,12] = gene_meta_fe$I2
    chr_meta_all_results[j,13] = r1_all_gene$BETA_Burden
    chr_meta_all_results[j,14] = r1_all_gene$SE_Burden
    chr_meta_all_results[j,15] = r1_all_gene$Pvalue_Burden
    chr_meta_all_results[j,16] = r2_all_gene$BETA_Burden
    chr_meta_all_results[j,17] = r2_all_gene$SE_Burden
    chr_meta_all_results[j,18] = r2_all_gene$Pvalue_Burden
  }
  all_meta = rbind(all_meta,chr_meta_all_results)
}

root_path <- "/lindstroem/austin_working/Dissertation/Meta/Aim1/cancer_groups"  
phenotype <- "ectoderm" 

meta_data_list <- list(
  pLoF_meta = "pLoF_meta",
  missense_meta = "missense_meta",
  synonymous_meta = "synonymous_meta",
  pLoF_missense_meta = "pLoF_missense_meta",
  missense_synonymous_meta = "missense_synonymous_meta",
  all_meta = "all_meta"
)

for (df_name in names(meta_data_list)) {
  file_path <- file.path(root_path, paste0(meta_data_list[[df_name]], "_", phenotype, ".txt"))
  fwrite(get(df_name), file_path, sep = "\t", quote = FALSE, na = "NA")
}

# Run mesoderm meta-analyses #

pLoF_meta = NULL
missense_meta = NULL
synonymous_meta = NULL
pLoF_missense_meta = NULL
missense_synonymous_meta = NULL
all_meta = NULL

for (i in 1:22) {
  UKBB_Ns = fread("/lindstroem/austin_working/Dissertation/UKBB/gene_based_results/cancer_groups/sparse_GRM/chr21_mesoderm_SAIGE_results.txt.singleAssoc.txt", nrows = 1) %>%
    select(N_case,N_ctrl)
  AoU_Ns = fread("/lindstroem/austin_working/Dissertation/AoU/gene_based_results/cancer_groups/sparse_GRM/chr21_mesoderm_SAIGE_results.txt.singleAssoc.txt", nrows = 1) %>%
    select(N_case,N_ctrl)
  Ns = data.frame(N_case = UKBB_Ns$N_case + AoU_Ns$N_case,
                  N_ctrl = UKBB_Ns$N_ctrl + AoU_Ns$N_ctrl) %>% 
    mutate(N_total = N_case + N_ctrl )
  
  
  UKBB_results = fread(paste0("/lindstroem/austin_working/Dissertation/UKBB/gene_based_results/cancer_groups/sparse_GRM/chr",i,"_mesoderm_SAIGE_results.txt"))
  AoU_results = fread(paste0("/lindstroem/austin_working/Dissertation/AoU/gene_based_results/cancer_groups/sparse_GRM/chr",i,"_mesoderm_SAIGE_results.txt"))
  
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
  
  
  chr_meta_pLoF_results = data.frame(Gene = NA,
                                     Group = NA,
                                     CHR = NA,
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
    
    chr_meta_pLoF_results[j,1] = common_genes_pLoF[j]
    chr_meta_pLoF_results[j,2] = r1_pLoF_gene$Group
    chr_meta_pLoF_results[j,3] = i
    chr_meta_pLoF_results[j,7] = gene_meta_fe$beta
    chr_meta_pLoF_results[j,8] = gene_meta_fe$se
    chr_meta_pLoF_results[j,9] = gene_meta_fe$ci.lb
    chr_meta_pLoF_results[j,10] = gene_meta_fe$ci.ub
    chr_meta_pLoF_results[j,11] = gene_meta_fe$pval
    chr_meta_pLoF_results[j,12] = gene_meta_fe$I2
    chr_meta_pLoF_results[j,13] = r1_pLoF_gene$BETA_Burden
    chr_meta_pLoF_results[j,14] = r1_pLoF_gene$SE_Burden
    chr_meta_pLoF_results[j,15] = r1_pLoF_gene$Pvalue_Burden
    chr_meta_pLoF_results[j,16] = r2_pLoF_gene$BETA_Burden
    chr_meta_pLoF_results[j,17] = r2_pLoF_gene$SE_Burden
    chr_meta_pLoF_results[j,18] = r2_pLoF_gene$Pvalue_Burden
  }
  pLoF_meta = rbind(pLoF_meta,chr_meta_pLoF_results)
  
  chr_meta_missense_results = data.frame(Gene = NA,
                                         Group = NA,
                                         CHR = NA,
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
    
    chr_meta_missense_results[j,1] = common_genes_missense[j]
    chr_meta_missense_results[j,2] = r1_missense_gene$Group
    chr_meta_missense_results[j,3] = i
    chr_meta_missense_results[j,7] = gene_meta_fe$beta
    chr_meta_missense_results[j,8] = gene_meta_fe$se
    chr_meta_missense_results[j,9] = gene_meta_fe$ci.lb
    chr_meta_missense_results[j,10] = gene_meta_fe$ci.ub
    chr_meta_missense_results[j,11] = gene_meta_fe$pval
    chr_meta_missense_results[j,12] = gene_meta_fe$I2
    chr_meta_missense_results[j,13] = r1_missense_gene$BETA_Burden
    chr_meta_missense_results[j,14] = r1_missense_gene$SE_Burden
    chr_meta_missense_results[j,15] = r1_missense_gene$Pvalue_Burden
    chr_meta_missense_results[j,16] = r2_missense_gene$BETA_Burden
    chr_meta_missense_results[j,17] = r2_missense_gene$SE_Burden
    chr_meta_missense_results[j,18] = r2_missense_gene$Pvalue_Burden
  }
  missense_meta = rbind(missense_meta,chr_meta_missense_results)
  
  chr_meta_pLoF_missense_results = data.frame(Gene = NA,
                                              Group = NA,
                                              CHR = NA,
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
    
    chr_meta_pLoF_missense_results[j,1] = common_genes_pLoF_missense[j]
    chr_meta_pLoF_missense_results[j,2] = r1_pLoF_missense_gene$Group
    chr_meta_pLoF_missense_results[j,3] = i
    chr_meta_pLoF_missense_results[j,7] = gene_meta_fe$beta
    chr_meta_pLoF_missense_results[j,8] = gene_meta_fe$se
    chr_meta_pLoF_missense_results[j,9] = gene_meta_fe$ci.lb
    chr_meta_pLoF_missense_results[j,10] = gene_meta_fe$ci.ub
    chr_meta_pLoF_missense_results[j,11] = gene_meta_fe$pval
    chr_meta_pLoF_missense_results[j,12] = gene_meta_fe$I2
    chr_meta_pLoF_missense_results[j,13] = r1_pLoF_missense_gene$BETA_Burden
    chr_meta_pLoF_missense_results[j,14] = r1_pLoF_missense_gene$SE_Burden
    chr_meta_pLoF_missense_results[j,15] = r1_pLoF_missense_gene$Pvalue_Burden
    chr_meta_pLoF_missense_results[j,16] = r2_pLoF_missense_gene$BETA_Burden
    chr_meta_pLoF_missense_results[j,17] = r2_pLoF_missense_gene$SE_Burden
    chr_meta_pLoF_missense_results[j,18] = r2_pLoF_missense_gene$Pvalue_Burden
  }
  pLoF_missense_meta = rbind(pLoF_missense_meta,chr_meta_pLoF_missense_results)
  
  chr_meta_missense_synonymous_results = data.frame(Gene = NA,
                                                    Group = NA,
                                                    CHR = NA,
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
    
    chr_meta_missense_synonymous_results[j,1] = common_genes_missense_synonymous[j]
    chr_meta_missense_synonymous_results[j,2] = r1_missense_synonymous_gene$Group
    chr_meta_missense_synonymous_results[j,3] = i
    chr_meta_missense_synonymous_results[j,7] = gene_meta_fe$beta
    chr_meta_missense_synonymous_results[j,8] = gene_meta_fe$se
    chr_meta_missense_synonymous_results[j,9] = gene_meta_fe$ci.lb
    chr_meta_missense_synonymous_results[j,10] = gene_meta_fe$ci.ub
    chr_meta_missense_synonymous_results[j,11] = gene_meta_fe$pval
    chr_meta_missense_synonymous_results[j,12] = gene_meta_fe$I2
    chr_meta_missense_synonymous_results[j,13] = r1_missense_synonymous_gene$BETA_Burden
    chr_meta_missense_synonymous_results[j,14] = r1_missense_synonymous_gene$SE_Burden
    chr_meta_missense_synonymous_results[j,15] = r1_missense_synonymous_gene$Pvalue_Burden
    chr_meta_missense_synonymous_results[j,16] = r2_missense_synonymous_gene$BETA_Burden
    chr_meta_missense_synonymous_results[j,17] = r2_missense_synonymous_gene$SE_Burden
    chr_meta_missense_synonymous_results[j,18] = r2_missense_synonymous_gene$Pvalue_Burden
  }
  missense_synonymous_meta = rbind(missense_synonymous_meta,chr_meta_missense_synonymous_results)
  
  chr_meta_synonymous_results = data.frame(Gene = NA,
                                           Group = NA,
                                           CHR = NA,
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
    
    chr_meta_synonymous_results[j,1] = common_genes_synonymous[j]
    chr_meta_synonymous_results[j,2] = r1_synonymous_gene$Group
    chr_meta_synonymous_results[j,3] = i
    chr_meta_synonymous_results[j,7] = gene_meta_fe$beta
    chr_meta_synonymous_results[j,8] = gene_meta_fe$se
    chr_meta_synonymous_results[j,9] = gene_meta_fe$ci.lb
    chr_meta_synonymous_results[j,10] = gene_meta_fe$ci.ub
    chr_meta_synonymous_results[j,11] = gene_meta_fe$pval
    chr_meta_synonymous_results[j,12] = gene_meta_fe$I2
    chr_meta_synonymous_results[j,13] = r1_synonymous_gene$BETA_Burden
    chr_meta_synonymous_results[j,14] = r1_synonymous_gene$SE_Burden
    chr_meta_synonymous_results[j,15] = r1_synonymous_gene$Pvalue_Burden
    chr_meta_synonymous_results[j,16] = r2_synonymous_gene$BETA_Burden
    chr_meta_synonymous_results[j,17] = r2_synonymous_gene$SE_Burden
    chr_meta_synonymous_results[j,18] = r2_synonymous_gene$Pvalue_Burden
  }
  synonymous_meta = rbind(synonymous_meta,chr_meta_synonymous_results)
  
  chr_meta_all_results = data.frame(Gene = NA,
                                    Group = NA,
                                    CHR = NA,
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
    
    chr_meta_all_results[j,1] = common_genes_all[j]
    chr_meta_all_results[j,2] = r1_all_gene$Group
    chr_meta_all_results[j,3] = i
    chr_meta_all_results[j,7] = gene_meta_fe$beta
    chr_meta_all_results[j,8] = gene_meta_fe$se
    chr_meta_all_results[j,9] = gene_meta_fe$ci.lb
    chr_meta_all_results[j,10] = gene_meta_fe$ci.ub
    chr_meta_all_results[j,11] = gene_meta_fe$pval
    chr_meta_all_results[j,12] = gene_meta_fe$I2
    chr_meta_all_results[j,13] = r1_all_gene$BETA_Burden
    chr_meta_all_results[j,14] = r1_all_gene$SE_Burden
    chr_meta_all_results[j,15] = r1_all_gene$Pvalue_Burden
    chr_meta_all_results[j,16] = r2_all_gene$BETA_Burden
    chr_meta_all_results[j,17] = r2_all_gene$SE_Burden
    chr_meta_all_results[j,18] = r2_all_gene$Pvalue_Burden
  }
  all_meta = rbind(all_meta,chr_meta_all_results)
}

root_path <- "/lindstroem/austin_working/Dissertation/Meta/Aim1/cancer_groups"  
phenotype <- "mesoderm" 

meta_data_list <- list(
  pLoF_meta = "pLoF_meta",
  missense_meta = "missense_meta",
  synonymous_meta = "synonymous_meta",
  pLoF_missense_meta = "pLoF_missense_meta",
  missense_synonymous_meta = "missense_synonymous_meta",
  all_meta = "all_meta"
)

for (df_name in names(meta_data_list)) {
  file_path <- file.path(root_path, paste0(meta_data_list[[df_name]], "_", phenotype, ".txt"))
  fwrite(get(df_name), file_path, sep = "\t", quote = FALSE, na = "NA")
}

# Run endoderm meta-analyses #

pLoF_meta = NULL
missense_meta = NULL
synonymous_meta = NULL
pLoF_missense_meta = NULL
missense_synonymous_meta = NULL
all_meta = NULL

for (i in 1:22) {
  UKBB_Ns = fread("/lindstroem/austin_working/Dissertation/UKBB/gene_based_results/cancer_groups/sparse_GRM/chr21_endoderm_SAIGE_results.txt.singleAssoc.txt", nrows = 1) %>%
    select(N_case,N_ctrl)
  AoU_Ns = fread("/lindstroem/austin_working/Dissertation/AoU/gene_based_results/cancer_groups/sparse_GRM/chr21_endoderm_SAIGE_results.txt.singleAssoc.txt", nrows = 1) %>%
    select(N_case,N_ctrl)
  Ns = data.frame(N_case = UKBB_Ns$N_case + AoU_Ns$N_case,
                  N_ctrl = UKBB_Ns$N_ctrl + AoU_Ns$N_ctrl) %>% 
    mutate(N_total = N_case + N_ctrl )
  
  
  UKBB_results = fread(paste0("/lindstroem/austin_working/Dissertation/UKBB/gene_based_results/cancer_groups/sparse_GRM/chr",i,"_endoderm_SAIGE_results.txt"))
  AoU_results = fread(paste0("/lindstroem/austin_working/Dissertation/AoU/gene_based_results/cancer_groups/sparse_GRM/chr",i,"_endoderm_SAIGE_results.txt"))
  
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
  
  
  chr_meta_pLoF_results = data.frame(Gene = NA,
                                     Group = NA,
                                     CHR = NA,
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
    
    chr_meta_pLoF_results[j,1] = common_genes_pLoF[j]
    chr_meta_pLoF_results[j,2] = r1_pLoF_gene$Group
    chr_meta_pLoF_results[j,3] = i
    chr_meta_pLoF_results[j,7] = gene_meta_fe$beta
    chr_meta_pLoF_results[j,8] = gene_meta_fe$se
    chr_meta_pLoF_results[j,9] = gene_meta_fe$ci.lb
    chr_meta_pLoF_results[j,10] = gene_meta_fe$ci.ub
    chr_meta_pLoF_results[j,11] = gene_meta_fe$pval
    chr_meta_pLoF_results[j,12] = gene_meta_fe$I2
    chr_meta_pLoF_results[j,13] = r1_pLoF_gene$BETA_Burden
    chr_meta_pLoF_results[j,14] = r1_pLoF_gene$SE_Burden
    chr_meta_pLoF_results[j,15] = r1_pLoF_gene$Pvalue_Burden
    chr_meta_pLoF_results[j,16] = r2_pLoF_gene$BETA_Burden
    chr_meta_pLoF_results[j,17] = r2_pLoF_gene$SE_Burden
    chr_meta_pLoF_results[j,18] = r2_pLoF_gene$Pvalue_Burden
  }
  pLoF_meta = rbind(pLoF_meta,chr_meta_pLoF_results)
  
  chr_meta_missense_results = data.frame(Gene = NA,
                                         Group = NA,
                                         CHR = NA,
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
    
    chr_meta_missense_results[j,1] = common_genes_missense[j]
    chr_meta_missense_results[j,2] = r1_missense_gene$Group
    chr_meta_missense_results[j,3] = i
    chr_meta_missense_results[j,7] = gene_meta_fe$beta
    chr_meta_missense_results[j,8] = gene_meta_fe$se
    chr_meta_missense_results[j,9] = gene_meta_fe$ci.lb
    chr_meta_missense_results[j,10] = gene_meta_fe$ci.ub
    chr_meta_missense_results[j,11] = gene_meta_fe$pval
    chr_meta_missense_results[j,12] = gene_meta_fe$I2
    chr_meta_missense_results[j,13] = r1_missense_gene$BETA_Burden
    chr_meta_missense_results[j,14] = r1_missense_gene$SE_Burden
    chr_meta_missense_results[j,15] = r1_missense_gene$Pvalue_Burden
    chr_meta_missense_results[j,16] = r2_missense_gene$BETA_Burden
    chr_meta_missense_results[j,17] = r2_missense_gene$SE_Burden
    chr_meta_missense_results[j,18] = r2_missense_gene$Pvalue_Burden
  }
  missense_meta = rbind(missense_meta,chr_meta_missense_results)
  
  chr_meta_pLoF_missense_results = data.frame(Gene = NA,
                                              Group = NA,
                                              CHR = NA,
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
    
    chr_meta_pLoF_missense_results[j,1] = common_genes_pLoF_missense[j]
    chr_meta_pLoF_missense_results[j,2] = r1_pLoF_missense_gene$Group
    chr_meta_pLoF_missense_results[j,3] = i
    chr_meta_pLoF_missense_results[j,7] = gene_meta_fe$beta
    chr_meta_pLoF_missense_results[j,8] = gene_meta_fe$se
    chr_meta_pLoF_missense_results[j,9] = gene_meta_fe$ci.lb
    chr_meta_pLoF_missense_results[j,10] = gene_meta_fe$ci.ub
    chr_meta_pLoF_missense_results[j,11] = gene_meta_fe$pval
    chr_meta_pLoF_missense_results[j,12] = gene_meta_fe$I2
    chr_meta_pLoF_missense_results[j,13] = r1_pLoF_missense_gene$BETA_Burden
    chr_meta_pLoF_missense_results[j,14] = r1_pLoF_missense_gene$SE_Burden
    chr_meta_pLoF_missense_results[j,15] = r1_pLoF_missense_gene$Pvalue_Burden
    chr_meta_pLoF_missense_results[j,16] = r2_pLoF_missense_gene$BETA_Burden
    chr_meta_pLoF_missense_results[j,17] = r2_pLoF_missense_gene$SE_Burden
    chr_meta_pLoF_missense_results[j,18] = r2_pLoF_missense_gene$Pvalue_Burden
  }
  pLoF_missense_meta = rbind(pLoF_missense_meta,chr_meta_pLoF_missense_results)
  
  chr_meta_missense_synonymous_results = data.frame(Gene = NA,
                                                    Group = NA,
                                                    CHR = NA,
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
    
    chr_meta_missense_synonymous_results[j,1] = common_genes_missense_synonymous[j]
    chr_meta_missense_synonymous_results[j,2] = r1_missense_synonymous_gene$Group
    chr_meta_missense_synonymous_results[j,3] = i
    chr_meta_missense_synonymous_results[j,7] = gene_meta_fe$beta
    chr_meta_missense_synonymous_results[j,8] = gene_meta_fe$se
    chr_meta_missense_synonymous_results[j,9] = gene_meta_fe$ci.lb
    chr_meta_missense_synonymous_results[j,10] = gene_meta_fe$ci.ub
    chr_meta_missense_synonymous_results[j,11] = gene_meta_fe$pval
    chr_meta_missense_synonymous_results[j,12] = gene_meta_fe$I2
    chr_meta_missense_synonymous_results[j,13] = r1_missense_synonymous_gene$BETA_Burden
    chr_meta_missense_synonymous_results[j,14] = r1_missense_synonymous_gene$SE_Burden
    chr_meta_missense_synonymous_results[j,15] = r1_missense_synonymous_gene$Pvalue_Burden
    chr_meta_missense_synonymous_results[j,16] = r2_missense_synonymous_gene$BETA_Burden
    chr_meta_missense_synonymous_results[j,17] = r2_missense_synonymous_gene$SE_Burden
    chr_meta_missense_synonymous_results[j,18] = r2_missense_synonymous_gene$Pvalue_Burden
  }
  missense_synonymous_meta = rbind(missense_synonymous_meta,chr_meta_missense_synonymous_results)
  
  chr_meta_synonymous_results = data.frame(Gene = NA,
                                           Group = NA,
                                           CHR = NA,
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
    
    chr_meta_synonymous_results[j,1] = common_genes_synonymous[j]
    chr_meta_synonymous_results[j,2] = r1_synonymous_gene$Group
    chr_meta_synonymous_results[j,3] = i
    chr_meta_synonymous_results[j,7] = gene_meta_fe$beta
    chr_meta_synonymous_results[j,8] = gene_meta_fe$se
    chr_meta_synonymous_results[j,9] = gene_meta_fe$ci.lb
    chr_meta_synonymous_results[j,10] = gene_meta_fe$ci.ub
    chr_meta_synonymous_results[j,11] = gene_meta_fe$pval
    chr_meta_synonymous_results[j,12] = gene_meta_fe$I2
    chr_meta_synonymous_results[j,13] = r1_synonymous_gene$BETA_Burden
    chr_meta_synonymous_results[j,14] = r1_synonymous_gene$SE_Burden
    chr_meta_synonymous_results[j,15] = r1_synonymous_gene$Pvalue_Burden
    chr_meta_synonymous_results[j,16] = r2_synonymous_gene$BETA_Burden
    chr_meta_synonymous_results[j,17] = r2_synonymous_gene$SE_Burden
    chr_meta_synonymous_results[j,18] = r2_synonymous_gene$Pvalue_Burden
  }
  synonymous_meta = rbind(synonymous_meta,chr_meta_synonymous_results)
  
  chr_meta_all_results = data.frame(Gene = NA,
                                    Group = NA,
                                    CHR = NA,
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
    
    chr_meta_all_results[j,1] = common_genes_all[j]
    chr_meta_all_results[j,2] = r1_all_gene$Group
    chr_meta_all_results[j,3] = i
    chr_meta_all_results[j,7] = gene_meta_fe$beta
    chr_meta_all_results[j,8] = gene_meta_fe$se
    chr_meta_all_results[j,9] = gene_meta_fe$ci.lb
    chr_meta_all_results[j,10] = gene_meta_fe$ci.ub
    chr_meta_all_results[j,11] = gene_meta_fe$pval
    chr_meta_all_results[j,12] = gene_meta_fe$I2
    chr_meta_all_results[j,13] = r1_all_gene$BETA_Burden
    chr_meta_all_results[j,14] = r1_all_gene$SE_Burden
    chr_meta_all_results[j,15] = r1_all_gene$Pvalue_Burden
    chr_meta_all_results[j,16] = r2_all_gene$BETA_Burden
    chr_meta_all_results[j,17] = r2_all_gene$SE_Burden
    chr_meta_all_results[j,18] = r2_all_gene$Pvalue_Burden
  }
  all_meta = rbind(all_meta,chr_meta_all_results)
}

root_path <- "/lindstroem/austin_working/Dissertation/Meta/Aim1/cancer_groups"  
phenotype <- "endoderm" 

meta_data_list <- list(
  pLoF_meta = "pLoF_meta",
  missense_meta = "missense_meta",
  synonymous_meta = "synonymous_meta",
  pLoF_missense_meta = "pLoF_missense_meta",
  missense_synonymous_meta = "missense_synonymous_meta",
  all_meta = "all_meta"
)

for (df_name in names(meta_data_list)) {
  file_path <- file.path(root_path, paste0(meta_data_list[[df_name]], "_", phenotype, ".txt"))
  fwrite(get(df_name), file_path, sep = "\t", quote = FALSE, na = "NA")
}

# Run hormone meta-analyses #

pLoF_meta = NULL
missense_meta = NULL
synonymous_meta = NULL
pLoF_missense_meta = NULL
missense_synonymous_meta = NULL
all_meta = NULL

for (i in 1:22) {
  UKBB_Ns = fread("/lindstroem/austin_working/Dissertation/UKBB/gene_based_results/cancer_groups/sparse_GRM/chr21_hormone_SAIGE_results.txt.singleAssoc.txt", nrows = 1) %>%
    select(N_case,N_ctrl)
  AoU_Ns = fread("/lindstroem/austin_working/Dissertation/AoU/gene_based_results/cancer_groups/sparse_GRM/chr21_hormone_SAIGE_results.txt.singleAssoc.txt", nrows = 1) %>%
    select(N_case,N_ctrl)
  Ns = data.frame(N_case = UKBB_Ns$N_case + AoU_Ns$N_case,
                  N_ctrl = UKBB_Ns$N_ctrl + AoU_Ns$N_ctrl) %>% 
    mutate(N_total = N_case + N_ctrl )
  
  
  UKBB_results = fread(paste0("/lindstroem/austin_working/Dissertation/UKBB/gene_based_results/cancer_groups/sparse_GRM/chr",i,"_hormone_SAIGE_results.txt"))
  AoU_results = fread(paste0("/lindstroem/austin_working/Dissertation/AoU/gene_based_results/cancer_groups/sparse_GRM/chr",i,"_hormone_SAIGE_results.txt"))
  
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
  
  
  chr_meta_pLoF_results = data.frame(Gene = NA,
                                     Group = NA,
                                     CHR = NA,
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
    
    chr_meta_pLoF_results[j,1] = common_genes_pLoF[j]
    chr_meta_pLoF_results[j,2] = r1_pLoF_gene$Group
    chr_meta_pLoF_results[j,3] = i
    chr_meta_pLoF_results[j,7] = gene_meta_fe$beta
    chr_meta_pLoF_results[j,8] = gene_meta_fe$se
    chr_meta_pLoF_results[j,9] = gene_meta_fe$ci.lb
    chr_meta_pLoF_results[j,10] = gene_meta_fe$ci.ub
    chr_meta_pLoF_results[j,11] = gene_meta_fe$pval
    chr_meta_pLoF_results[j,12] = gene_meta_fe$I2
    chr_meta_pLoF_results[j,13] = r1_pLoF_gene$BETA_Burden
    chr_meta_pLoF_results[j,14] = r1_pLoF_gene$SE_Burden
    chr_meta_pLoF_results[j,15] = r1_pLoF_gene$Pvalue_Burden
    chr_meta_pLoF_results[j,16] = r2_pLoF_gene$BETA_Burden
    chr_meta_pLoF_results[j,17] = r2_pLoF_gene$SE_Burden
    chr_meta_pLoF_results[j,18] = r2_pLoF_gene$Pvalue_Burden
  }
  pLoF_meta = rbind(pLoF_meta,chr_meta_pLoF_results)
  
  chr_meta_missense_results = data.frame(Gene = NA,
                                         Group = NA,
                                         CHR = NA,
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
    
    chr_meta_missense_results[j,1] = common_genes_missense[j]
    chr_meta_missense_results[j,2] = r1_missense_gene$Group
    chr_meta_missense_results[j,3] = i
    chr_meta_missense_results[j,7] = gene_meta_fe$beta
    chr_meta_missense_results[j,8] = gene_meta_fe$se
    chr_meta_missense_results[j,9] = gene_meta_fe$ci.lb
    chr_meta_missense_results[j,10] = gene_meta_fe$ci.ub
    chr_meta_missense_results[j,11] = gene_meta_fe$pval
    chr_meta_missense_results[j,12] = gene_meta_fe$I2
    chr_meta_missense_results[j,13] = r1_missense_gene$BETA_Burden
    chr_meta_missense_results[j,14] = r1_missense_gene$SE_Burden
    chr_meta_missense_results[j,15] = r1_missense_gene$Pvalue_Burden
    chr_meta_missense_results[j,16] = r2_missense_gene$BETA_Burden
    chr_meta_missense_results[j,17] = r2_missense_gene$SE_Burden
    chr_meta_missense_results[j,18] = r2_missense_gene$Pvalue_Burden
  }
  missense_meta = rbind(missense_meta,chr_meta_missense_results)
  
  chr_meta_pLoF_missense_results = data.frame(Gene = NA,
                                              Group = NA,
                                              CHR = NA,
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
    
    chr_meta_pLoF_missense_results[j,1] = common_genes_pLoF_missense[j]
    chr_meta_pLoF_missense_results[j,2] = r1_pLoF_missense_gene$Group
    chr_meta_pLoF_missense_results[j,3] = i
    chr_meta_pLoF_missense_results[j,7] = gene_meta_fe$beta
    chr_meta_pLoF_missense_results[j,8] = gene_meta_fe$se
    chr_meta_pLoF_missense_results[j,9] = gene_meta_fe$ci.lb
    chr_meta_pLoF_missense_results[j,10] = gene_meta_fe$ci.ub
    chr_meta_pLoF_missense_results[j,11] = gene_meta_fe$pval
    chr_meta_pLoF_missense_results[j,12] = gene_meta_fe$I2
    chr_meta_pLoF_missense_results[j,13] = r1_pLoF_missense_gene$BETA_Burden
    chr_meta_pLoF_missense_results[j,14] = r1_pLoF_missense_gene$SE_Burden
    chr_meta_pLoF_missense_results[j,15] = r1_pLoF_missense_gene$Pvalue_Burden
    chr_meta_pLoF_missense_results[j,16] = r2_pLoF_missense_gene$BETA_Burden
    chr_meta_pLoF_missense_results[j,17] = r2_pLoF_missense_gene$SE_Burden
    chr_meta_pLoF_missense_results[j,18] = r2_pLoF_missense_gene$Pvalue_Burden
  }
  pLoF_missense_meta = rbind(pLoF_missense_meta,chr_meta_pLoF_missense_results)
  
  chr_meta_missense_synonymous_results = data.frame(Gene = NA,
                                                    Group = NA,
                                                    CHR = NA,
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
    
    chr_meta_missense_synonymous_results[j,1] = common_genes_missense_synonymous[j]
    chr_meta_missense_synonymous_results[j,2] = r1_missense_synonymous_gene$Group
    chr_meta_missense_synonymous_results[j,3] = i
    chr_meta_missense_synonymous_results[j,7] = gene_meta_fe$beta
    chr_meta_missense_synonymous_results[j,8] = gene_meta_fe$se
    chr_meta_missense_synonymous_results[j,9] = gene_meta_fe$ci.lb
    chr_meta_missense_synonymous_results[j,10] = gene_meta_fe$ci.ub
    chr_meta_missense_synonymous_results[j,11] = gene_meta_fe$pval
    chr_meta_missense_synonymous_results[j,12] = gene_meta_fe$I2
    chr_meta_missense_synonymous_results[j,13] = r1_missense_synonymous_gene$BETA_Burden
    chr_meta_missense_synonymous_results[j,14] = r1_missense_synonymous_gene$SE_Burden
    chr_meta_missense_synonymous_results[j,15] = r1_missense_synonymous_gene$Pvalue_Burden
    chr_meta_missense_synonymous_results[j,16] = r2_missense_synonymous_gene$BETA_Burden
    chr_meta_missense_synonymous_results[j,17] = r2_missense_synonymous_gene$SE_Burden
    chr_meta_missense_synonymous_results[j,18] = r2_missense_synonymous_gene$Pvalue_Burden
  }
  missense_synonymous_meta = rbind(missense_synonymous_meta,chr_meta_missense_synonymous_results)
  
  chr_meta_synonymous_results = data.frame(Gene = NA,
                                           Group = NA,
                                           CHR = NA,
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
    
    chr_meta_synonymous_results[j,1] = common_genes_synonymous[j]
    chr_meta_synonymous_results[j,2] = r1_synonymous_gene$Group
    chr_meta_synonymous_results[j,3] = i
    chr_meta_synonymous_results[j,7] = gene_meta_fe$beta
    chr_meta_synonymous_results[j,8] = gene_meta_fe$se
    chr_meta_synonymous_results[j,9] = gene_meta_fe$ci.lb
    chr_meta_synonymous_results[j,10] = gene_meta_fe$ci.ub
    chr_meta_synonymous_results[j,11] = gene_meta_fe$pval
    chr_meta_synonymous_results[j,12] = gene_meta_fe$I2
    chr_meta_synonymous_results[j,13] = r1_synonymous_gene$BETA_Burden
    chr_meta_synonymous_results[j,14] = r1_synonymous_gene$SE_Burden
    chr_meta_synonymous_results[j,15] = r1_synonymous_gene$Pvalue_Burden
    chr_meta_synonymous_results[j,16] = r2_synonymous_gene$BETA_Burden
    chr_meta_synonymous_results[j,17] = r2_synonymous_gene$SE_Burden
    chr_meta_synonymous_results[j,18] = r2_synonymous_gene$Pvalue_Burden
  }
  synonymous_meta = rbind(synonymous_meta,chr_meta_synonymous_results)
  
  chr_meta_all_results = data.frame(Gene = NA,
                                    Group = NA,
                                    CHR = NA,
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
    
    chr_meta_all_results[j,1] = common_genes_all[j]
    chr_meta_all_results[j,2] = r1_all_gene$Group
    chr_meta_all_results[j,3] = i
    chr_meta_all_results[j,7] = gene_meta_fe$beta
    chr_meta_all_results[j,8] = gene_meta_fe$se
    chr_meta_all_results[j,9] = gene_meta_fe$ci.lb
    chr_meta_all_results[j,10] = gene_meta_fe$ci.ub
    chr_meta_all_results[j,11] = gene_meta_fe$pval
    chr_meta_all_results[j,12] = gene_meta_fe$I2
    chr_meta_all_results[j,13] = r1_all_gene$BETA_Burden
    chr_meta_all_results[j,14] = r1_all_gene$SE_Burden
    chr_meta_all_results[j,15] = r1_all_gene$Pvalue_Burden
    chr_meta_all_results[j,16] = r2_all_gene$BETA_Burden
    chr_meta_all_results[j,17] = r2_all_gene$SE_Burden
    chr_meta_all_results[j,18] = r2_all_gene$Pvalue_Burden
  }
  all_meta = rbind(all_meta,chr_meta_all_results)
}

root_path <- "/lindstroem/austin_working/Dissertation/Meta/Aim1/cancer_groups"  
phenotype <- "hormone" 

meta_data_list <- list(
  pLoF_meta = "pLoF_meta",
  missense_meta = "missense_meta",
  synonymous_meta = "synonymous_meta",
  pLoF_missense_meta = "pLoF_missense_meta",
  missense_synonymous_meta = "missense_synonymous_meta",
  all_meta = "all_meta"
)

for (df_name in names(meta_data_list)) {
  file_path <- file.path(root_path, paste0(meta_data_list[[df_name]], "_", phenotype, ".txt"))
  fwrite(get(df_name), file_path, sep = "\t", quote = FALSE, na = "NA")
}

# Run infectious meta-analyses #

pLoF_meta = NULL
missense_meta = NULL
synonymous_meta = NULL
pLoF_missense_meta = NULL
missense_synonymous_meta = NULL
all_meta = NULL

for (i in 1:22) {
  UKBB_Ns = fread("/lindstroem/austin_working/Dissertation/UKBB/gene_based_results/cancer_groups/sparse_GRM/chr21_infectious_SAIGE_results.txt.singleAssoc.txt", nrows = 1) %>%
    select(N_case,N_ctrl)
  AoU_Ns = fread("/lindstroem/austin_working/Dissertation/AoU/gene_based_results/cancer_groups/sparse_GRM/chr21_infectious_SAIGE_results.txt.singleAssoc.txt", nrows = 1) %>%
    select(N_case,N_ctrl)
  Ns = data.frame(N_case = UKBB_Ns$N_case + AoU_Ns$N_case,
                  N_ctrl = UKBB_Ns$N_ctrl + AoU_Ns$N_ctrl) %>% 
    mutate(N_total = N_case + N_ctrl )
  
  
  UKBB_results = fread(paste0("/lindstroem/austin_working/Dissertation/UKBB/gene_based_results/cancer_groups/sparse_GRM/chr",i,"_infectious_SAIGE_results.txt"))
  AoU_results = fread(paste0("/lindstroem/austin_working/Dissertation/AoU/gene_based_results/cancer_groups/sparse_GRM/chr",i,"_infectious_SAIGE_results.txt"))
  
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
  
  
  chr_meta_pLoF_results = data.frame(Gene = NA,
                                     Group = NA,
                                     CHR = NA,
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
    
    chr_meta_pLoF_results[j,1] = common_genes_pLoF[j]
    chr_meta_pLoF_results[j,2] = r1_pLoF_gene$Group
    chr_meta_pLoF_results[j,3] = i
    chr_meta_pLoF_results[j,7] = gene_meta_fe$beta
    chr_meta_pLoF_results[j,8] = gene_meta_fe$se
    chr_meta_pLoF_results[j,9] = gene_meta_fe$ci.lb
    chr_meta_pLoF_results[j,10] = gene_meta_fe$ci.ub
    chr_meta_pLoF_results[j,11] = gene_meta_fe$pval
    chr_meta_pLoF_results[j,12] = gene_meta_fe$I2
    chr_meta_pLoF_results[j,13] = r1_pLoF_gene$BETA_Burden
    chr_meta_pLoF_results[j,14] = r1_pLoF_gene$SE_Burden
    chr_meta_pLoF_results[j,15] = r1_pLoF_gene$Pvalue_Burden
    chr_meta_pLoF_results[j,16] = r2_pLoF_gene$BETA_Burden
    chr_meta_pLoF_results[j,17] = r2_pLoF_gene$SE_Burden
    chr_meta_pLoF_results[j,18] = r2_pLoF_gene$Pvalue_Burden
  }
  pLoF_meta = rbind(pLoF_meta,chr_meta_pLoF_results)
  
  chr_meta_missense_results = data.frame(Gene = NA,
                                         Group = NA,
                                         CHR = NA,
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
    
    chr_meta_missense_results[j,1] = common_genes_missense[j]
    chr_meta_missense_results[j,2] = r1_missense_gene$Group
    chr_meta_missense_results[j,3] = i
    chr_meta_missense_results[j,7] = gene_meta_fe$beta
    chr_meta_missense_results[j,8] = gene_meta_fe$se
    chr_meta_missense_results[j,9] = gene_meta_fe$ci.lb
    chr_meta_missense_results[j,10] = gene_meta_fe$ci.ub
    chr_meta_missense_results[j,11] = gene_meta_fe$pval
    chr_meta_missense_results[j,12] = gene_meta_fe$I2
    chr_meta_missense_results[j,13] = r1_missense_gene$BETA_Burden
    chr_meta_missense_results[j,14] = r1_missense_gene$SE_Burden
    chr_meta_missense_results[j,15] = r1_missense_gene$Pvalue_Burden
    chr_meta_missense_results[j,16] = r2_missense_gene$BETA_Burden
    chr_meta_missense_results[j,17] = r2_missense_gene$SE_Burden
    chr_meta_missense_results[j,18] = r2_missense_gene$Pvalue_Burden
  }
  missense_meta = rbind(missense_meta,chr_meta_missense_results)
  
  chr_meta_pLoF_missense_results = data.frame(Gene = NA,
                                              Group = NA,
                                              CHR = NA,
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
    
    chr_meta_pLoF_missense_results[j,1] = common_genes_pLoF_missense[j]
    chr_meta_pLoF_missense_results[j,2] = r1_pLoF_missense_gene$Group
    chr_meta_pLoF_missense_results[j,3] = i
    chr_meta_pLoF_missense_results[j,7] = gene_meta_fe$beta
    chr_meta_pLoF_missense_results[j,8] = gene_meta_fe$se
    chr_meta_pLoF_missense_results[j,9] = gene_meta_fe$ci.lb
    chr_meta_pLoF_missense_results[j,10] = gene_meta_fe$ci.ub
    chr_meta_pLoF_missense_results[j,11] = gene_meta_fe$pval
    chr_meta_pLoF_missense_results[j,12] = gene_meta_fe$I2
    chr_meta_pLoF_missense_results[j,13] = r1_pLoF_missense_gene$BETA_Burden
    chr_meta_pLoF_missense_results[j,14] = r1_pLoF_missense_gene$SE_Burden
    chr_meta_pLoF_missense_results[j,15] = r1_pLoF_missense_gene$Pvalue_Burden
    chr_meta_pLoF_missense_results[j,16] = r2_pLoF_missense_gene$BETA_Burden
    chr_meta_pLoF_missense_results[j,17] = r2_pLoF_missense_gene$SE_Burden
    chr_meta_pLoF_missense_results[j,18] = r2_pLoF_missense_gene$Pvalue_Burden
  }
  pLoF_missense_meta = rbind(pLoF_missense_meta,chr_meta_pLoF_missense_results)
  
  chr_meta_missense_synonymous_results = data.frame(Gene = NA,
                                                    Group = NA,
                                                    CHR = NA,
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
    
    chr_meta_missense_synonymous_results[j,1] = common_genes_missense_synonymous[j]
    chr_meta_missense_synonymous_results[j,2] = r1_missense_synonymous_gene$Group
    chr_meta_missense_synonymous_results[j,3] = i
    chr_meta_missense_synonymous_results[j,7] = gene_meta_fe$beta
    chr_meta_missense_synonymous_results[j,8] = gene_meta_fe$se
    chr_meta_missense_synonymous_results[j,9] = gene_meta_fe$ci.lb
    chr_meta_missense_synonymous_results[j,10] = gene_meta_fe$ci.ub
    chr_meta_missense_synonymous_results[j,11] = gene_meta_fe$pval
    chr_meta_missense_synonymous_results[j,12] = gene_meta_fe$I2
    chr_meta_missense_synonymous_results[j,13] = r1_missense_synonymous_gene$BETA_Burden
    chr_meta_missense_synonymous_results[j,14] = r1_missense_synonymous_gene$SE_Burden
    chr_meta_missense_synonymous_results[j,15] = r1_missense_synonymous_gene$Pvalue_Burden
    chr_meta_missense_synonymous_results[j,16] = r2_missense_synonymous_gene$BETA_Burden
    chr_meta_missense_synonymous_results[j,17] = r2_missense_synonymous_gene$SE_Burden
    chr_meta_missense_synonymous_results[j,18] = r2_missense_synonymous_gene$Pvalue_Burden
  }
  missense_synonymous_meta = rbind(missense_synonymous_meta,chr_meta_missense_synonymous_results)
  
  chr_meta_synonymous_results = data.frame(Gene = NA,
                                           Group = NA,
                                           CHR = NA,
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
    
    chr_meta_synonymous_results[j,1] = common_genes_synonymous[j]
    chr_meta_synonymous_results[j,2] = r1_synonymous_gene$Group
    chr_meta_synonymous_results[j,3] = i
    chr_meta_synonymous_results[j,7] = gene_meta_fe$beta
    chr_meta_synonymous_results[j,8] = gene_meta_fe$se
    chr_meta_synonymous_results[j,9] = gene_meta_fe$ci.lb
    chr_meta_synonymous_results[j,10] = gene_meta_fe$ci.ub
    chr_meta_synonymous_results[j,11] = gene_meta_fe$pval
    chr_meta_synonymous_results[j,12] = gene_meta_fe$I2
    chr_meta_synonymous_results[j,13] = r1_synonymous_gene$BETA_Burden
    chr_meta_synonymous_results[j,14] = r1_synonymous_gene$SE_Burden
    chr_meta_synonymous_results[j,15] = r1_synonymous_gene$Pvalue_Burden
    chr_meta_synonymous_results[j,16] = r2_synonymous_gene$BETA_Burden
    chr_meta_synonymous_results[j,17] = r2_synonymous_gene$SE_Burden
    chr_meta_synonymous_results[j,18] = r2_synonymous_gene$Pvalue_Burden
  }
  synonymous_meta = rbind(synonymous_meta,chr_meta_synonymous_results)
  
  chr_meta_all_results = data.frame(Gene = NA,
                                    Group = NA,
                                    CHR = NA,
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
    
    chr_meta_all_results[j,1] = common_genes_all[j]
    chr_meta_all_results[j,2] = r1_all_gene$Group
    chr_meta_all_results[j,3] = i
    chr_meta_all_results[j,7] = gene_meta_fe$beta
    chr_meta_all_results[j,8] = gene_meta_fe$se
    chr_meta_all_results[j,9] = gene_meta_fe$ci.lb
    chr_meta_all_results[j,10] = gene_meta_fe$ci.ub
    chr_meta_all_results[j,11] = gene_meta_fe$pval
    chr_meta_all_results[j,12] = gene_meta_fe$I2
    chr_meta_all_results[j,13] = r1_all_gene$BETA_Burden
    chr_meta_all_results[j,14] = r1_all_gene$SE_Burden
    chr_meta_all_results[j,15] = r1_all_gene$Pvalue_Burden
    chr_meta_all_results[j,16] = r2_all_gene$BETA_Burden
    chr_meta_all_results[j,17] = r2_all_gene$SE_Burden
    chr_meta_all_results[j,18] = r2_all_gene$Pvalue_Burden
  }
  all_meta = rbind(all_meta,chr_meta_all_results)
}

root_path <- "/lindstroem/austin_working/Dissertation/Meta/Aim1/cancer_groups"  
phenotype <- "infectious" 

meta_data_list <- list(
  pLoF_meta = "pLoF_meta",
  missense_meta = "missense_meta",
  synonymous_meta = "synonymous_meta",
  pLoF_missense_meta = "pLoF_missense_meta",
  missense_synonymous_meta = "missense_synonymous_meta",
  all_meta = "all_meta"
)

for (df_name in names(meta_data_list)) {
  file_path <- file.path(root_path, paste0(meta_data_list[[df_name]], "_", phenotype, ".txt"))
  fwrite(get(df_name), file_path, sep = "\t", quote = FALSE, na = "NA")
}

# Run smoking meta-analyses #

pLoF_meta = NULL
missense_meta = NULL
synonymous_meta = NULL
pLoF_missense_meta = NULL
missense_synonymous_meta = NULL
all_meta = NULL

for (i in 1:22) {
  UKBB_Ns = fread("/lindstroem/austin_working/Dissertation/UKBB/gene_based_results/cancer_groups/sparse_GRM/chr21_smoking_SAIGE_results.txt.singleAssoc.txt", nrows = 1) %>%
    select(N_case,N_ctrl)
  AoU_Ns = fread("/lindstroem/austin_working/Dissertation/AoU/gene_based_results/cancer_groups/sparse_GRM/chr21_smoking_SAIGE_results.txt.singleAssoc.txt", nrows = 1) %>%
    select(N_case,N_ctrl)
  Ns = data.frame(N_case = UKBB_Ns$N_case + AoU_Ns$N_case,
                  N_ctrl = UKBB_Ns$N_ctrl + AoU_Ns$N_ctrl) %>% 
    mutate(N_total = N_case + N_ctrl )
  
  
  UKBB_results = fread(paste0("/lindstroem/austin_working/Dissertation/UKBB/gene_based_results/cancer_groups/sparse_GRM/chr",i,"_smoking_SAIGE_results.txt"))
  AoU_results = fread(paste0("/lindstroem/austin_working/Dissertation/AoU/gene_based_results/cancer_groups/sparse_GRM/chr",i,"_smoking_SAIGE_results.txt"))
  
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
  
  
  chr_meta_pLoF_results = data.frame(Gene = NA,
                                     Group = NA,
                                     CHR = NA,
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
    
    chr_meta_pLoF_results[j,1] = common_genes_pLoF[j]
    chr_meta_pLoF_results[j,2] = r1_pLoF_gene$Group
    chr_meta_pLoF_results[j,3] = i
    chr_meta_pLoF_results[j,7] = gene_meta_fe$beta
    chr_meta_pLoF_results[j,8] = gene_meta_fe$se
    chr_meta_pLoF_results[j,9] = gene_meta_fe$ci.lb
    chr_meta_pLoF_results[j,10] = gene_meta_fe$ci.ub
    chr_meta_pLoF_results[j,11] = gene_meta_fe$pval
    chr_meta_pLoF_results[j,12] = gene_meta_fe$I2
    chr_meta_pLoF_results[j,13] = r1_pLoF_gene$BETA_Burden
    chr_meta_pLoF_results[j,14] = r1_pLoF_gene$SE_Burden
    chr_meta_pLoF_results[j,15] = r1_pLoF_gene$Pvalue_Burden
    chr_meta_pLoF_results[j,16] = r2_pLoF_gene$BETA_Burden
    chr_meta_pLoF_results[j,17] = r2_pLoF_gene$SE_Burden
    chr_meta_pLoF_results[j,18] = r2_pLoF_gene$Pvalue_Burden
  }
  pLoF_meta = rbind(pLoF_meta,chr_meta_pLoF_results)
  
  chr_meta_missense_results = data.frame(Gene = NA,
                                         Group = NA,
                                         CHR = NA,
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
    
    chr_meta_missense_results[j,1] = common_genes_missense[j]
    chr_meta_missense_results[j,2] = r1_missense_gene$Group
    chr_meta_missense_results[j,3] = i
    chr_meta_missense_results[j,7] = gene_meta_fe$beta
    chr_meta_missense_results[j,8] = gene_meta_fe$se
    chr_meta_missense_results[j,9] = gene_meta_fe$ci.lb
    chr_meta_missense_results[j,10] = gene_meta_fe$ci.ub
    chr_meta_missense_results[j,11] = gene_meta_fe$pval
    chr_meta_missense_results[j,12] = gene_meta_fe$I2
    chr_meta_missense_results[j,13] = r1_missense_gene$BETA_Burden
    chr_meta_missense_results[j,14] = r1_missense_gene$SE_Burden
    chr_meta_missense_results[j,15] = r1_missense_gene$Pvalue_Burden
    chr_meta_missense_results[j,16] = r2_missense_gene$BETA_Burden
    chr_meta_missense_results[j,17] = r2_missense_gene$SE_Burden
    chr_meta_missense_results[j,18] = r2_missense_gene$Pvalue_Burden
  }
  missense_meta = rbind(missense_meta,chr_meta_missense_results)
  
  chr_meta_pLoF_missense_results = data.frame(Gene = NA,
                                              Group = NA,
                                              CHR = NA,
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
    
    chr_meta_pLoF_missense_results[j,1] = common_genes_pLoF_missense[j]
    chr_meta_pLoF_missense_results[j,2] = r1_pLoF_missense_gene$Group
    chr_meta_pLoF_missense_results[j,3] = i
    chr_meta_pLoF_missense_results[j,7] = gene_meta_fe$beta
    chr_meta_pLoF_missense_results[j,8] = gene_meta_fe$se
    chr_meta_pLoF_missense_results[j,9] = gene_meta_fe$ci.lb
    chr_meta_pLoF_missense_results[j,10] = gene_meta_fe$ci.ub
    chr_meta_pLoF_missense_results[j,11] = gene_meta_fe$pval
    chr_meta_pLoF_missense_results[j,12] = gene_meta_fe$I2
    chr_meta_pLoF_missense_results[j,13] = r1_pLoF_missense_gene$BETA_Burden
    chr_meta_pLoF_missense_results[j,14] = r1_pLoF_missense_gene$SE_Burden
    chr_meta_pLoF_missense_results[j,15] = r1_pLoF_missense_gene$Pvalue_Burden
    chr_meta_pLoF_missense_results[j,16] = r2_pLoF_missense_gene$BETA_Burden
    chr_meta_pLoF_missense_results[j,17] = r2_pLoF_missense_gene$SE_Burden
    chr_meta_pLoF_missense_results[j,18] = r2_pLoF_missense_gene$Pvalue_Burden
  }
  pLoF_missense_meta = rbind(pLoF_missense_meta,chr_meta_pLoF_missense_results)
  
  chr_meta_missense_synonymous_results = data.frame(Gene = NA,
                                                    Group = NA,
                                                    CHR = NA,
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
    
    chr_meta_missense_synonymous_results[j,1] = common_genes_missense_synonymous[j]
    chr_meta_missense_synonymous_results[j,2] = r1_missense_synonymous_gene$Group
    chr_meta_missense_synonymous_results[j,3] = i
    chr_meta_missense_synonymous_results[j,7] = gene_meta_fe$beta
    chr_meta_missense_synonymous_results[j,8] = gene_meta_fe$se
    chr_meta_missense_synonymous_results[j,9] = gene_meta_fe$ci.lb
    chr_meta_missense_synonymous_results[j,10] = gene_meta_fe$ci.ub
    chr_meta_missense_synonymous_results[j,11] = gene_meta_fe$pval
    chr_meta_missense_synonymous_results[j,12] = gene_meta_fe$I2
    chr_meta_missense_synonymous_results[j,13] = r1_missense_synonymous_gene$BETA_Burden
    chr_meta_missense_synonymous_results[j,14] = r1_missense_synonymous_gene$SE_Burden
    chr_meta_missense_synonymous_results[j,15] = r1_missense_synonymous_gene$Pvalue_Burden
    chr_meta_missense_synonymous_results[j,16] = r2_missense_synonymous_gene$BETA_Burden
    chr_meta_missense_synonymous_results[j,17] = r2_missense_synonymous_gene$SE_Burden
    chr_meta_missense_synonymous_results[j,18] = r2_missense_synonymous_gene$Pvalue_Burden
  }
  missense_synonymous_meta = rbind(missense_synonymous_meta,chr_meta_missense_synonymous_results)
  
  chr_meta_synonymous_results = data.frame(Gene = NA,
                                           Group = NA,
                                           CHR = NA,
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
    
    chr_meta_synonymous_results[j,1] = common_genes_synonymous[j]
    chr_meta_synonymous_results[j,2] = r1_synonymous_gene$Group
    chr_meta_synonymous_results[j,3] = i
    chr_meta_synonymous_results[j,7] = gene_meta_fe$beta
    chr_meta_synonymous_results[j,8] = gene_meta_fe$se
    chr_meta_synonymous_results[j,9] = gene_meta_fe$ci.lb
    chr_meta_synonymous_results[j,10] = gene_meta_fe$ci.ub
    chr_meta_synonymous_results[j,11] = gene_meta_fe$pval
    chr_meta_synonymous_results[j,12] = gene_meta_fe$I2
    chr_meta_synonymous_results[j,13] = r1_synonymous_gene$BETA_Burden
    chr_meta_synonymous_results[j,14] = r1_synonymous_gene$SE_Burden
    chr_meta_synonymous_results[j,15] = r1_synonymous_gene$Pvalue_Burden
    chr_meta_synonymous_results[j,16] = r2_synonymous_gene$BETA_Burden
    chr_meta_synonymous_results[j,17] = r2_synonymous_gene$SE_Burden
    chr_meta_synonymous_results[j,18] = r2_synonymous_gene$Pvalue_Burden
  }
  synonymous_meta = rbind(synonymous_meta,chr_meta_synonymous_results)
  
  chr_meta_all_results = data.frame(Gene = NA,
                                    Group = NA,
                                    CHR = NA,
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
    
    chr_meta_all_results[j,1] = common_genes_all[j]
    chr_meta_all_results[j,2] = r1_all_gene$Group
    chr_meta_all_results[j,3] = i
    chr_meta_all_results[j,7] = gene_meta_fe$beta
    chr_meta_all_results[j,8] = gene_meta_fe$se
    chr_meta_all_results[j,9] = gene_meta_fe$ci.lb
    chr_meta_all_results[j,10] = gene_meta_fe$ci.ub
    chr_meta_all_results[j,11] = gene_meta_fe$pval
    chr_meta_all_results[j,12] = gene_meta_fe$I2
    chr_meta_all_results[j,13] = r1_all_gene$BETA_Burden
    chr_meta_all_results[j,14] = r1_all_gene$SE_Burden
    chr_meta_all_results[j,15] = r1_all_gene$Pvalue_Burden
    chr_meta_all_results[j,16] = r2_all_gene$BETA_Burden
    chr_meta_all_results[j,17] = r2_all_gene$SE_Burden
    chr_meta_all_results[j,18] = r2_all_gene$Pvalue_Burden
  }
  all_meta = rbind(all_meta,chr_meta_all_results)
}

root_path <- "/lindstroem/austin_working/Dissertation/Meta/Aim1/cancer_groups"  
phenotype <- "smoking" 

meta_data_list <- list(
  pLoF_meta = "pLoF_meta",
  missense_meta = "missense_meta",
  synonymous_meta = "synonymous_meta",
  pLoF_missense_meta = "pLoF_missense_meta",
  missense_synonymous_meta = "missense_synonymous_meta",
  all_meta = "all_meta"
)

for (df_name in names(meta_data_list)) {
  file_path <- file.path(root_path, paste0(meta_data_list[[df_name]], "_", phenotype, ".txt"))
  fwrite(get(df_name), file_path, sep = "\t", quote = FALSE, na = "NA")
}



