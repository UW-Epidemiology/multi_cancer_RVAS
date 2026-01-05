# Setup #
rm(list = ls())
install.packages("pacman")
library(pacman)
p_load("tidyverse","data.table")


# Load and prepare cleaned phenotype data #

## UKB cleaned core phenotype data ##
system("dx download file-GyfJxgQJfyQQYFJ9f874Qbpj")
UKB_core = fread("UKB_core_phenotypes_2_2025_cleaned_dataset.csv",
                 check.names = TRUE, na.strings = "")

## UKB cleaned cancer phenotype data #
system("dx download file-GyfPkF8JfyQVkfJ0ZB91v79P")
UKB_cancer = fread("UKB_cancer_phenotypes_2_2025_cleaned_dataset.csv",
                   check.names = TRUE, na.strings = "")

## Load participant relatedness data ##
system("dx download file-Gjfx320Jqkj0BffVXjQ92QFp")
UKB_rel = fread("ukb_rel.dat",
                check.names = TRUE, na.strings = "")


## Merge core and cancer phenotype datasets ##
UKB_cancer = merge(UKB_cancer,
                   UKB_core, 
                   by="Participant_ID",
                   all.x = TRUE)
rm(UKB_core)

# Exclude individuals with discordant self-reported and genetically inferred sex #
UKB_cancer = UKB_cancer %>%
  filter(Sex == Genetic_sex)

# Filter the dataset to remove pairs with inferred 1st degree or closer relationships #

## Set the relatedness threshold at 1st degree or closer ##
rel_thresh = (1/(2^(5/2)))
table(UKB_rel$Kinship > rel_thresh)

## Extract a list of case and control ids ##
case_ids = UKB_cancer %>%
  filter(CancerCase==1) %>%
  pull(Participant_ID) %>%
  unique()

control_ids = UKB_cancer %>%
  filter(CancerCase==0) %>%
  pull(Participant_ID) %>%
  unique()


## Define optimized function to remove individuals with high relatedness (preferentially dropping controls) ##
remove_related_individuals <- function(df, cases, controls, rel_thresh) {
  df <- as.data.frame(df)
  # Filter out rows where kinship coefficient exceeds the threshold
  high_kinship_pairs <- df[df[, 5] > rel_thresh, ]
  id1 <- high_kinship_pairs[, 1]
  id2 <- high_kinship_pairs[, 2]
  id1_is_case <- id1 %in% cases
  id2_is_case <- id2 %in% cases
  id1_is_control <- id1 %in% controls
  id2_is_control <- id2 %in% controls
  removed_ids <- vector()
  
  # Remove controls preferentially, otherwise remove one arbitrarily
  removed_ids <- c(removed_ids, 
                   id1[id1_is_control & !id2_is_control],
                   id2[id2_is_control & !id1_is_control],
                   id1[id1_is_control & id2_is_control],
                   id1[id1_is_case & id2_is_case],
                   id2[id1_is_case & id2_is_control],
                   id1[id2_is_case & id1_is_control])
  
  # Remove duplicates
  removed_ids <- unique(removed_ids)
  
  # Determine the final list of IDs to keep
  ids_to_keep <- setdiff(c(cases, controls), removed_ids)
  
  # Report the number of cases and controls before and after exclusions
  final_cases <- intersect(ids_to_keep, cases)
  final_controls <- intersect(ids_to_keep, controls)
  
  cat("Original number of cases:", length(cases), "\n")
  cat("Original number of controls:", length(controls), "\n")
  cat("Number of cases removed:", length(cases) - length(final_cases), "\n")
  cat("Number of controls removed:", length(controls) - length(final_controls), "\n")
  cat("Number of cases after exclusions:", length(final_cases), "\n")
  cat("Number of controls after exclusions:", length(final_controls), "\n")
  
  # Return the vector of participant identifiers to keep
  return(ids_to_keep)
}

## Drop one sample from each pair with estimated kinship of 1st degree or closer ##
set.seed(134689)
ids_to_keep = remove_related_individuals(UKB_rel,case_ids,control_ids,rel_thresh)
cancer_data = UKB_cancer %>%
  filter(Participant_ID %in% ids_to_keep)

# # Write files with participant IDs for analyses #
# 
# ## IDs (single column per row) ##
# IDs = cancer_data %>% 
#   pull(Participant_ID) %>%
#   unique()
# 
# ## PLINK-format IDs (two columns per row) ##
# PLINK_format_IDs = data.frame(FID = IDs,IID = IDs)
# 
# 
# ## Write and upload ID files ##
# write.table(IDs, file = "samples.txt",
#             row.names = FALSE,
#             col.names = FALSE,
#             sep = " ",
#             quote = FALSE)
# write.table(PLINK_format_IDs,
#             file = "PLINK_samples.txt",
#             row.names = FALSE,
#             col.names = FALSE,
#             sep = " ",
#             quote = FALSE)

# # Save the sample filtering files #
# system("dx upload samples.txt --path project-GZqFkYQJfyQzB1G6YBg89ZG3:/austin_working/dissertation_project/cancer_phenotypes/")
# system("dx upload PLINK_samples.txt --path project-GZqFkYQJfyQzB1G6YBg89ZG3:/austin_working/dissertation_project/cancer_phenotypes/")

# Filter to IDs in WES dataset #
system("dx download file-Gjfx5g8JqkjB3GKvq7zFBK79")
WES_IDs = fread("ukb23158_c1_b0_v1.fam")
cancer_data = cancer_data %>%
  filter(Participant_ID %in% WES_IDs$V1)

# Define case and control sets #

## Case group definitions ##

ectoderm = c("lip", "tongue", "gum", "mouth floor", "palate",
             "other mouth", "parotid gland", "salivary glands", "accessory sinuses",
             "thymus", "melanoma", "non-melanoma skin", "peripheral nervous system", "breast",
             "prostate", "eye", "meninges", "brain", "other central nervous")
mesoderm = c("heart and other", "bone and cartilage", "mesothelioma", "Kaposi's sarcoma",
             "retroperitoneum and peritoneum", "other soft tissue", "cervix", "other uterus",
             "endometrium", "ovary", "placenta", "testis", "kidney", "renal pelvis", "ureter",
             "bladder", "Hodgkin's lymphoma", "non-follicular lymphoma",
             "mature T and NK-cell lymphoma", "other non-Hodgkin's lymphoma",
             "other T and NK-cell lymphoma", "malignant immunoproliferative disease",
             "multiple myeloma", "lymphoid leukemia", "myeloid leukemia", "monocytic leukemia",
             "other leukemia")
endoderm = c("tonsil", "oropharynx", "nasopharynx", "pyriform sinus", "hypopharynx", "esophagus",
             "stomach", "small intestine", "appendix", "colon", "rectosigmoid junction", "rectum",
             "anus", "liver", "gallbladder", "pancreas", "larynx", "trachea", "lung and bronchus",
             "vagina", "thyroid gland")
hormone = c("breast","endometrium","ovary","prostate","testis","thyroid gland")
smoking = c("lip", "tongue", "gum", "mouth floor", "palate", "other mouth", "parotid gland",
            "salivary glands", "tonsil", "oropharynx", "nasopharynx", "pyriform sinus", "hypopharynx",
            "other head and neck", "esophagus", "stomach", "liver", "pancreas", "larynx", "trachea",
            "lung and bronchus", "cervix", "kidney", "renal pelvis", "ureter", "bladder",
            "colon", "rectosigmoid junction", "rectum","anus")
infectious = c("lip", "tongue", "gum", "mouth floor", "palate", "other mouth", "parotid gland",
               "salivary glands", "tonsil", "oropharynx", "nasopharynx", "pyriform sinus", "hypopharynx",
               "other head and neck", "stomach", "anus", "liver", "gallbladder", "larynx",
               "lung and bronchus", "Kaposi's sarcoma", "vulva", "cervix", "Hodgkin's lymphoma",
               "follicular lymphoma", "non-follicular lymphoma", "mature T and NK-cell lymphoma",
               "other non-Hodgkin's lymphoma","other T and NK-cell lymphoma")

## Add a female genetic sex indicator variable ##
cancer_data = cancer_data %>%
  mutate(Female = ifelse(Genetic_sex == "Female",1,0))

## Subset cancer cases by cancer type or diagnosis group##
in_ectoderm = cancer_data %>%
  filter(CancerSite %in% ectoderm) %>%
  select(Participant_ID)

in_mesoderm = cancer_data %>%
  filter(CancerSite %in% mesoderm) %>%
  select(Participant_ID)

in_endoderm = cancer_data %>%
  filter(CancerSite %in% endoderm) %>%
  select(Participant_ID)

in_hormone = cancer_data %>%
  filter(CancerSite %in% hormone) %>%
  select(Participant_ID)

in_smoking = cancer_data %>%
  filter(CancerSite %in% smoking) %>%
  select(Participant_ID)

in_infectious = cancer_data %>%
  filter(CancerSite %in% infectious) %>%
  select(Participant_ID)

in_anus = cancer_data %>%
  filter(CancerSite == "anus") %>%
  select(Participant_ID)

in_bladder = cancer_data %>%
  filter(CancerSite == "bladder") %>%
  select(Participant_ID)

in_bone = cancer_data %>%
  filter(CancerSite == "bone and cartilage") %>%
  select(Participant_ID)

in_brain = cancer_data %>%
  filter(CancerSite == "brain") %>%
  select(Participant_ID)

in_breast = cancer_data %>%
  filter(CancerSite == "breast") %>%
  select(Participant_ID)

in_cervix = cancer_data %>%
  filter(CancerSite == "cervix", Female == 1) %>%
  select(Participant_ID)

in_colon = cancer_data %>%
  filter(CancerSite == "colon") %>%
  select(Participant_ID)

in_endometrium = cancer_data %>%
  filter(CancerSite == "endometrium", Female == 1) %>%
  select(Participant_ID)

in_esophagus = cancer_data %>%
  filter(CancerSite == "esophagus") %>%
  select(Participant_ID)

in_eye = cancer_data %>%
  filter(CancerSite == "eye") %>%
  select(Participant_ID)

in_follicular_lymphoma = cancer_data %>%
  filter(CancerSite == "follicular lymphoma") %>%
  select(Participant_ID)

in_heart = cancer_data %>%
  filter(CancerSite == "heart and other") %>%
  select(Participant_ID)

in_Hodgkins_lymphoma = cancer_data %>%
  filter(CancerSite == "Hodgkin's lymphoma") %>%
  select(Participant_ID)

in_kidney = cancer_data %>%
  filter(CancerSite == "kidney") %>%
  select(Participant_ID)

in_larynx = cancer_data %>%
  filter(CancerSite == "larynx") %>%
  select(Participant_ID)

in_liver = cancer_data %>%
  filter(CancerSite == "liver") %>%
  select(Participant_ID)

in_lung = cancer_data %>%
  filter(CancerSite == "lung and bronchus") %>%
  select(Participant_ID)

in_lymphoid_leukemia = cancer_data %>%
  filter(CancerSite == "lymphoid leukemia") %>%
  select(Participant_ID)

in_MID = cancer_data %>%
  filter(CancerSite == "malignant immunoproliferative disease") %>%
  select(Participant_ID)

in_M_TNK_lymphoma = cancer_data %>%
  filter(CancerSite == "mature T and NK-cell lymphoma") %>%
  select(Participant_ID)

in_melanoma = cancer_data %>%
  filter(CancerSite == "melanoma") %>%
  select(Participant_ID)

in_myeloma = cancer_data %>%
  filter(CancerSite == "multiple myeloma") %>%
  select(Participant_ID)

in_myeloid_leukemia = cancer_data %>%
  filter(CancerSite == "myeloid leukemia") %>%
  select(Participant_ID)

in_NF_lymphoma = cancer_data %>%
  filter(CancerSite == "non-follicular lymphoma") %>%
  select(Participant_ID)

in_non_melanoma = cancer_data %>%
  filter(CancerSite == "non-melanoma skin") %>%
  select(Participant_ID)

in_oropharynx = cancer_data %>%
  filter(CancerSite == "oropharynx") %>%
  select(Participant_ID)

in_OCNS = cancer_data %>%
  filter(CancerSite == "other soft tissue") %>%
  select(Participant_ID)

in_other_uterus = cancer_data %>%
  filter(CancerSite == "other uterus", Female == 1) %>%
  select(Participant_ID)

in_ovary = cancer_data %>%
  filter(CancerSite == "ovary", Female == 1) %>%
  select(Participant_ID)

in_pancreas = cancer_data %>%
  filter(CancerSite == "pancreas") %>%
  select(Participant_ID)

in_prostate = cancer_data %>%
  filter(CancerSite == "prostate", Female == 0) %>%
  select(Participant_ID)

in_rectosigmoid = cancer_data %>%
  filter(CancerSite == "rectosigmoid junction") %>%
  select(Participant_ID)

in_rectum = cancer_data %>%
  filter(CancerSite == "rectum") %>%
  select(Participant_ID)

in_renal_pelvis = cancer_data %>%
  filter(CancerSite == "renal pelvis") %>%
  select(Participant_ID)

in_peritoneum = cancer_data %>%
  filter(CancerSite == "retroperitoneum and peritoneum") %>%
  select(Participant_ID)

in_small_intestine = cancer_data %>%
  filter(CancerSite == "small intestine") %>%
  select(Participant_ID)

in_stomach = cancer_data %>%
  filter(CancerSite == "stomach") %>%
  select(Participant_ID)

in_testis = cancer_data %>%
  filter(CancerSite == "testis", Female == 0) %>%
  select(Participant_ID)

in_thyroid = cancer_data %>%
  filter(CancerSite == "thyroid gland") %>%
  select(Participant_ID)

in_tongue = cancer_data %>%
  filter(CancerSite == "tongue") %>%
  select(Participant_ID)

in_tonsil = cancer_data %>%
  filter(CancerSite == "tonsil") %>%
  select(Participant_ID)

## Test case-control sets ##

# ### Cases ##
# case_data = cancer_data %>% filter(CancerCase == 1)
# 
# test_breast_pheno = case_data %>%
#   filter(Participant_ID %in% in_breast$Participant_ID & Female == 1) %>% 
#   filter(CancerSite %in% c("breast")) %>%
#   group_by(Participant_ID) %>%
#   slice_min(DxAge) %>%
#   ungroup %>%
#   select(c("Participant_ID",
#            "Age_at_recruitment",
#            paste0("UKB_PC",1:20))) %>%
#   distinct()
# 
# test_endometrium_pheno = case_data %>%
#   filter(Participant_ID %in% in_endometrium$Participant_ID & Female == 1) %>% 
#   filter(CancerSite %in% c("endometrium")) %>%
#   group_by(Participant_ID) %>%
#   slice_min(DxAge) %>%
#   ungroup %>%
#   select(c("Participant_ID",
#            "Age_at_recruitment",
#            paste0("UKB_PC",1:20))) %>%
#   distinct()
# 
# test_ovary_pheno = case_data %>%
#   filter(Participant_ID %in% in_ovary$Participant_ID & Female == 1) %>% 
#   filter(CancerSite %in% c("ovary")) %>%
#   group_by(Participant_ID) %>%
#   slice_min(DxAge) %>%
#   ungroup %>%
#   select(c("Participant_ID",
#            "Age_at_recruitment",
#            paste0("UKB_PC",1:20))) %>%
#   distinct()
# 
# ### Resolve case sample overlaps ###
# resolve_sample_overlap <- function(dataset1, dataset2, dataset3, id_col = "Participant_ID") {
#   # Get dataset sizes
#   sizes <- c(nrow(dataset1), nrow(dataset2), nrow(dataset3))
#   
#   # Extract unique sample IDs
#   ids <- list(dataset1[[id_col]], dataset2[[id_col]], dataset3[[id_col]])
#   
#   # Identify overlapping sample IDs
#   common_12 <- intersect(ids[[1]], ids[[2]])
#   common_13 <- intersect(ids[[1]], ids[[3]])
#   common_23 <- intersect(ids[[2]], ids[[3]])
#   common_all <- Reduce(intersect, ids)
#   
#   # Order the datasets by size (largest first)
#   dataset_order <- order(sizes, decreasing = TRUE)
#   
#   # Create a list of the datasets to modify by size order
#   datasets <- list(dataset1, dataset2, dataset3)
#   datasets <- datasets[dataset_order]
#   
#   # Apply filtering in size order, ensuring larger datasets lose the overlaps first
#   datasets[[1]] <- datasets[[1]] %>% filter(!Participant_ID %in% common_all & !Participant_ID %in% common_12 & !Participant_ID %in% common_13)
#   datasets[[2]] <- datasets[[2]] %>% filter(!Participant_ID %in% common_all & !Participant_ID %in% common_23)
#   datasets[[3]] <- datasets[[3]]
#   
#   # Return the updated datasets
#   dataset1 <- datasets[[1]]
#   dataset2 <- datasets[[2]]
#   dataset3 <- datasets[[3]]
#   
#   # Direct assignment of updated datasets to original names
#   assign("dataset1", dataset1, envir = parent.frame())
#   assign("dataset2", dataset2, envir = parent.frame())
#   assign("dataset3", dataset3, envir = parent.frame())
# }
# 
# 
# resolve_sample_overlap(test_breast_pheno,
#                        test_endometrium_pheno,
#                        test_ovary_pheno)
# 
# test_breast_pheno = dataset1
# test_endometrium_pheno = dataset2
# test_ovary_pheno = dataset3
# 
# ### Controls ###
# control_data = cancer_data %>% filter(CancerCase == 0)
# 
# test_control_pheno = control_data %>%
#   select(c("Participant_ID",
#            "Age_at_recruitment",
#            paste0("UKB_PC",1:20)))
# 
# set.seed(216156)
# num_breast_control_rows <- nrow(test_breast_pheno) * 2 
# num_endometrium_control_rows <- nrow(test_endometrium_pheno) * 5  
# num_ovary_control_rows <- nrow(test_ovary_pheno) * 5  
# 
# test_control_breast_pheno = test_control_pheno %>%
#   slice_sample(n = min(num_breast_control_rows, nrow(test_control_pheno)))
# test_control_endometrium_pheno = test_control_pheno %>%
#   slice_sample(n = min(num_endometrium_control_rows, nrow(test_control_pheno)))
# test_control_ovary_pheno = test_control_pheno %>%
#   slice_sample(n = min(num_ovary_control_rows, nrow(test_control_pheno)))
# 
# 
# 
# ### Resolve control sample overlaps ###
# resolve_sample_overlap(test_control_breast_pheno,
#                        test_control_endometrium_pheno,
#                        test_control_ovary_pheno)
# 
# test_control_breast_pheno = dataset1
# test_control_endometrium_pheno = dataset2
# test_control_ovary_pheno = dataset3
# 
# ### Merge test case and control sets ##
# test_breast_pheno = test_breast_pheno %>% mutate(breast = 1, endometrium = NA, ovary = NA)
# test_endometrium_pheno = test_endometrium_pheno %>% mutate(breast = NA, endometrium = 1, ovary = NA)
# test_ovary_pheno = test_ovary_pheno %>% mutate(breast = NA, endometrium =NA, ovary = 1)
# test_control_breast_pheno = test_control_breast_pheno %>% mutate(breast = 0, endometrium = NA, ovary = NA)
# test_control_endometrium_pheno = test_control_endometrium_pheno %>% mutate(breast = NA, endometrium = 0, ovary = NA)
# test_control_ovary_pheno = test_control_ovary_pheno %>% mutate(breast = NA, endometrium = NA, ovary = 0)
# 
# test_pheno = rbind(test_breast_pheno,
#                    test_endometrium_pheno,
#                    test_ovary_pheno,
#                    test_control_breast_pheno,
#                    test_control_endometrium_pheno,
#                    test_control_ovary_pheno)
# 
# ## Write test phenotype data ##
# write.table(test_pheno,
#             file = "test_phenotypes.txt",
#             row.names = FALSE,
#             col.names = TRUE,
#             sep = " ",
#             quote = FALSE)
# 
# # Save the sample filtering files #
# system("dx upload test_phenotypes.txt --path project-GZqFkYQJfyQzB1G6YBg89ZG3:/austin_working/dissertation_project/cancer_phenotypes/")

## Case-control sets ##
case_data = cancer_data %>%
  filter(CancerCase == 1)
unique_cases <- !duplicated(case_data$Participant_ID)
unique_cases <- case_data[unique_cases,]

case_pheno = unique_cases %>%
  select(c("Participant_ID",
           "Age_at_recruitment",
           "Female",
            paste0("UKB_PC",1:20)))

case_pheno = case_pheno %>%
  mutate(ectoderm = ifelse(Participant_ID %in% in_ectoderm$Participant_ID,1,NA),
         mesoderm = ifelse(Participant_ID %in% in_mesoderm$Participant_ID,1,NA),
         endoderm = ifelse(Participant_ID %in% in_endoderm$Participant_ID,1,NA),
         hormone = ifelse(Participant_ID %in% in_hormone$Participant_ID,1,NA),
         smoking = ifelse(Participant_ID %in% in_smoking$Participant_ID,1,NA),
         infectious = ifelse(Participant_ID %in% in_infectious$Participant_ID,1,NA),
         anus = ifelse(Participant_ID %in% in_anus$Participant_ID,1,NA),
         bladder = ifelse(Participant_ID %in% in_bladder$Participant_ID,1,NA),
         brain = ifelse(Participant_ID %in% in_brain$Participant_ID,1,NA),
         breast = ifelse(Participant_ID %in% in_breast$Participant_ID,1,NA),
         cervix = ifelse(Participant_ID %in% in_cervix$Participant_ID,1,NA),
         colon = ifelse(Participant_ID %in% in_colon$Participant_ID,1,NA),
         endometrium = ifelse(Participant_ID %in% in_endometrium$Participant_ID,1,NA),
         esophagus = ifelse(Participant_ID %in% in_esophagus$Participant_ID,1,NA),
         eye = ifelse(Participant_ID %in% in_eye$Participant_ID,1,NA),
         follicular_lymphoma = ifelse(Participant_ID %in% in_follicular_lymphoma$Participant_ID,1,NA),
         heart = ifelse(Participant_ID %in% in_heart$Participant_ID,1,NA),
         HL = ifelse(Participant_ID %in% in_Hodgkins_lymphoma$Participant_ID,1,NA),
         kidney = ifelse(Participant_ID %in% in_kidney$Participant_ID,1,NA),
         larynx = ifelse(Participant_ID %in% in_larynx$Participant_ID,1,NA),
         liver = ifelse(Participant_ID %in% in_liver$Participant_ID,1,NA),
         lung = ifelse(Participant_ID %in% in_lung$Participant_ID,1,NA),
         lymphoid_leukemia = ifelse(Participant_ID %in% in_lymphoid_leukemia$Participant_ID,1,NA),
         MID = ifelse(Participant_ID %in% in_MID$Participant_ID,1,NA),
         M_TNK_lymphoma = ifelse(Participant_ID %in% in_M_TNK_lymphoma$Participant_ID,1,NA),
         myeloid_leukemia = ifelse(Participant_ID %in% in_myeloid_leukemia$Participant_ID,1,NA),
         NF_lymphoma = ifelse(Participant_ID %in% in_NF_lymphoma$Participant_ID,1,NA),
         non_melanoma = ifelse(Participant_ID %in% in_non_melanoma$Participant_ID,1,NA),
         oropharynx = ifelse(Participant_ID %in% in_oropharynx$Participant_ID,1,NA),
         OCNS = ifelse(Participant_ID %in% in_OCNS$Participant_ID,1,NA),
         other_uterus = ifelse(Participant_ID %in% in_other_uterus$Participant_ID,1,NA),
         ovary = ifelse(Participant_ID %in% in_ovary$Participant_ID,1,NA),
         pancreas = ifelse(Participant_ID %in% in_pancreas$Participant_ID,1,NA),
         prostate = ifelse(Participant_ID %in% in_prostate$Participant_ID,1,NA),
         rectosigmoid = ifelse(Participant_ID %in% in_rectosigmoid$Participant_ID,1,NA),
         rectum = ifelse(Participant_ID %in% in_rectum$Participant_ID,1,NA),
         renal_pelvis = ifelse(Participant_ID %in% in_renal_pelvis$Participant_ID,1,NA),
         peritoneum = ifelse(Participant_ID %in% in_peritoneum$Participant_ID,1,NA),
         small_intestine = ifelse(Participant_ID %in% in_small_intestine$Participant_ID,1,NA),
         stomach = ifelse(Participant_ID %in% in_stomach$Participant_ID,1,NA),
         testis = ifelse(Participant_ID %in% in_testis$Participant_ID,1,NA),
         thyroid = ifelse(Participant_ID %in% in_thyroid$Participant_ID,1,NA),
         tongue = ifelse(Participant_ID %in% in_tongue$Participant_ID,1,NA),
         tonsil = ifelse(Participant_ID %in% in_tonsil$Participant_ID,1,NA))

control_data = cancer_data %>%
  filter(CancerCase == 0)

control_pheno = control_data %>%
  select(c("Participant_ID",
           "Age_at_recruitment",
           "Female",
           paste0("UKB_PC",1:20)))

control_pheno = control_pheno %>%
  mutate(ectoderm = 0,
         mesoderm = 0,
         endoderm = 0,
         hormone = 0,
         smoking = 0,
         infectious = 0,
         anus = 0,
         bladder = 0,
         brain = 0,
         breast = ifelse(Female == 0, NA, 0),
         cervix = ifelse(Female == 0, NA, 0),
         colon = 0,
         endometrium = ifelse(Female == 0, NA, 0),
         esophagus = 0,
         eye = 0,
         follicular_lymphoma = 0,
         heart = 0,
         HL = 0,
         kidney = 0,
         larynx = 0,
         liver = 0,
         lung = 0,
         lymphoid_leukemia = 0,
         MID = 0,
         M_TNK_lymphoma = 0,
         myeloid_leukemia = 0,
         NF_lymphoma = 0,
         non_melanoma = 0,
         oropharynx = 0,
         OCNS = 0,
         other_uterus = ifelse(Female == 0, NA, 0),
         ovary = ifelse(Female == 0, NA, 0),
         pancreas = 0,
         prostate = ifelse(Female == 0, 0, NA),
         rectosigmoid = 0,
         rectum = 0,
         renal_pelvis = 0,
         peritoneum = 0,
         small_intestine = 0,
         stomach = 0,
         testis = ifelse(Female == 0, 0, NA),
         thyroid = 0,
         tongue = 0,
         tonsil = 0)

## Combine into final phenotype dataset ##
pheno = rbind(case_pheno,control_pheno) %>%
  arrange(Participant_ID)

## Write final phenotype dataset ##
write.table(pheno,
            file = "phenotypes.txt",
            row.names = FALSE,
            col.names = TRUE,
            sep = " ",
            quote = FALSE)

## Update ##

# Save dataset and R script #
system("dx upload phenotypes.txt --path project-GZqFkYQJfyQzB1G6YBg89ZG3:/austin_working/dissertation_project/cancer_phenotypes/")

system("dx upload define_cancer_group_phenotype_files.R --path project-GZqFkYQJfyQzB1G6YBg89ZG3:/austin_working/dissertation_project/cancer_phenotypes/scripts/")

