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

# Filter to IDs in WES dataset #
system("dx download file-Gjfx5g8JqkjB3GKvq7zFBK79")
WES_IDs = fread("ukb23158_c1_b0_v1.fam")
cancer_data = cancer_data %>%
  filter(Participant_ID %in% WES_IDs$V1)

# Define case and control sets #

## Small case group definitions ##

colorectal = c("colon","rectum","rectosigmoid junction")
NHL = c("non-follicular lymphoma","mature T and NK-cell lymphoma", "other non-Hodgkin's lymphoma",
             "other T and NK-cell lymphoma", "malignant immunoproliferative disease",
             "multiple myeloma")
leukemias = c("lymphoid leukemia", "myeloid leukemia", "monocytic leukemia",
             "other leukemia")
oral = c("lip", "tongue", "gum", "mouth floor", "palate",
             "other mouth", "parotid gland", "salivary glands", "accessory sinuses")
neck = c("oropharynx", "nasopharynx", "pyriform sinus", "hypopharynx","larynx")

## Add a female genetic sex indicator variable ##
cancer_data = cancer_data %>%
  mutate(Female = ifelse(Genetic_sex == "Female",1,0))

## Subset cancer cases by cancer type or small group##
in_colorectal = cancer_data %>%
  filter(CancerSite %in% colorectal) %>%
  select(Participant_ID)

in_NHL = cancer_data %>%
  filter(CancerSite %in% NHL) %>%
  select(Participant_ID)

in_leukemias = cancer_data %>%
  filter(CancerSite %in% leukemias) %>%
  select(Participant_ID)

in_oral = cancer_data %>%
  filter(CancerSite %in% oral) %>%
  select(Participant_ID)

in_neck = cancer_data %>%
  filter(CancerSite %in% neck) %>%
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
  filter(CancerSite == "breast", Female == 1) %>%
  select(Participant_ID)

in_cervix = cancer_data %>%
  filter(CancerSite == "cervix", Female == 1) %>%
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

in_Hodgkins_lymphoma = cancer_data %>%
  filter(CancerSite == "Hodgkin's lymphoma") %>%
  select(Participant_ID)

in_kidney = cancer_data %>%
  filter(CancerSite == "kidney") %>%
  select(Participant_ID)

in_liver = cancer_data %>%
  filter(CancerSite == "liver") %>%
  select(Participant_ID)

in_lung = cancer_data %>%
  filter(CancerSite == "lung and bronchus") %>%
  select(Participant_ID)

in_melanoma = cancer_data %>%
  filter(CancerSite == "melanoma") %>%
  select(Participant_ID)

in_myeloma = cancer_data %>%
  filter(CancerSite == "multiple myeloma") %>%
  select(Participant_ID)

in_non_melanoma = cancer_data %>%
  filter(CancerSite == "non-melanoma skin") %>%
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

in_stomach = cancer_data %>%
  filter(CancerSite == "stomach") %>%
  select(Participant_ID)

in_testis = cancer_data %>%
  filter(CancerSite == "testis", Female == 0) %>%
  select(Participant_ID)

in_thyroid = cancer_data %>%
  filter(CancerSite == "thyroid gland") %>%
  select(Participant_ID)

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
  mutate(colorectal = ifelse(Participant_ID %in% in_colorectal$Participant_ID,1,NA),
         NHL = ifelse(Participant_ID %in% in_NHL$Participant_ID,1,NA),
         leukemias = ifelse(Participant_ID %in% in_leukemias$Participant_ID,1,NA),
         oral = ifelse(Participant_ID %in% in_oral$Participant_ID,1,NA),
         neck = ifelse(Participant_ID %in% in_neck$Participant_ID,1,NA),
         anus = ifelse(Participant_ID %in% in_anus$Participant_ID,1,NA),
         bladder = ifelse(Participant_ID %in% in_bladder$Participant_ID,1,NA),
         brain = ifelse(Participant_ID %in% in_brain$Participant_ID,1,NA),
         breast = ifelse(Participant_ID %in% in_breast$Participant_ID,1,NA),
         bone = ifelse(Participant_ID %in% in_bone$Participant_ID,1,NA),
         prostate = ifelse(Participant_ID %in% in_prostate$Participant_ID,1,NA),
         cervix = ifelse(Participant_ID %in% in_cervix$Participant_ID,1,NA),
         endometrium = ifelse(Participant_ID %in% in_endometrium$Participant_ID,1,NA),
         esophagus = ifelse(Participant_ID %in% in_esophagus$Participant_ID,1,NA),
         eye = ifelse(Participant_ID %in% in_eye$Participant_ID,1,NA),
         HL = ifelse(Participant_ID %in% in_Hodgkins_lymphoma$Participant_ID,1,NA),
         kidney = ifelse(Participant_ID %in% in_kidney$Participant_ID,1,NA),
         liver = ifelse(Participant_ID %in% in_liver$Participant_ID,1,NA),
         lung = ifelse(Participant_ID %in% in_lung$Participant_ID,1,NA),
         non_melanoma = ifelse(Participant_ID %in% in_non_melanoma$Participant_ID,1,NA),
         melanoma = ifelse(Participant_ID %in% in_melanoma$Participant_ID,1,NA),
         myeloma = ifelse(Participant_ID %in% in_myeloma$Participant_ID,1,NA),
         ovary = ifelse(Participant_ID %in% in_ovary$Participant_ID,1,NA),
         pancreas = ifelse(Participant_ID %in% in_pancreas$Participant_ID,1,NA),
         stomach = ifelse(Participant_ID %in% in_stomach$Participant_ID,1,NA),
         testis = ifelse(Participant_ID %in% in_testis$Participant_ID,1,NA),
         thyroid = ifelse(Participant_ID %in% in_thyroid$Participant_ID,1,NA))

control_data = cancer_data %>%
  filter(CancerCase == 0)

control_pheno = control_data %>%
  select(c("Participant_ID",
           "Age_at_recruitment",
           "Female",
           paste0("UKB_PC",1:20)))

control_pheno = control_pheno %>%
  mutate(colorectal = 0,
         NHL = 0,
         leukemias = 0,
         oral = 0,
         neck = 0,
         anus = 0,
         bladder = 0,
         brain = 0,
         breast = ifelse(Female == 0, NA, 0),
         bone = 0,
         prostate = ifelse(Female == 0, 0, NA),
         cervix = ifelse(Female == 0, NA, 0),
         endometrium = ifelse(Female == 0, NA, 0),
         esophagus = 0,
         eye = 0,
         HL = 0,
         kidney = 0,
         liver = 0,
         lung = 0,
         non_melanoma = 0,
         melanoma = 0,
         myeloma = 0,
         ovary = ifelse(Female == 0, NA, 0),
         pancreas = 0,
         stomach = 0,
         testis = ifelse(Female == 0, 0, NA),
         thyroid = 0)

## Combine into final phenotype dataset ##
pheno = rbind(case_pheno,control_pheno) %>%
  arrange(Participant_ID)

## Write final phenotype dataset ##
write.table(pheno,
            file = "single_cancer_phenotypes.txt",
            row.names = FALSE,
            col.names = TRUE,
            sep = " ",
            quote = FALSE)

## Update ##

# Save dataset and R script #
system("dx upload single_cancer_phenotypes.txt --path project-GZqFkYQJfyQzB1G6YBg89ZG3:/austin_working/dissertation_project/cancer_phenotypes/")

system("dx upload define_individual_cancer_phenotype_files.R --path project-GZqFkYQJfyQzB1G6YBg89ZG3:/austin_working/dissertation_project/cancer_phenotypes/scripts/")

