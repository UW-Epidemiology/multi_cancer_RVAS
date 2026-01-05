# Setup #
rm(list = ls())
library(tidyverse)
library(data.table)

# Prepare loop variables #
cancers = c("anus", "bladder", "bone", "breast", "brain", "cervix",
            "colorectal", "endometrium", "esophagus", "eye", "HL",
            "kidney", "leukemias", "liver", "lung", "melanoma",
            "myeloma", "NHL", "neck", "non_melanoma", "oral", "ovary",
            "pancreas", "prostate", "stomach", "testis", "thyroid")

# # Loop to consolidate results for each cancer analysis #
# for (i in cancers) {
#       chr_gene_results = fread(paste0("/lindstroem/austin_working/Dissertation/AoU/gene_based_results/individual_cancers_combined/",i,"_SAIGE_results.txt"))
#       chr_single_results = fread(paste0("/lindstroem/austin_working/Dissertation/AoU/gene_based_results/individual_cancers_combined/",i,"_SAIGE_results.txt.singleAssoc.txt"))
#       SAMHD1_gene_results = fread(paste0("/lindstroem/austin_working/Dissertation/AoU/gene_based_results/individual_cancers_combined/SAMHD1/chr20_",i,"_SAIGE_results.txt"))
#       SAMHD1_single_results = fread(paste0("/lindstroem/austin_working/Dissertation/AoU/gene_based_results/individual_cancers_combined/SAMHD1/chr20_",i,"_SAIGE_results.txt.singleAssoc.txt"))
#       gene_results = rbind.data.frame(chr_gene_results,SAMHD1_gene_results)
#       single_results = rbind.data.frame( chr_single_results,SAMHD1_single_results)
#       fwrite(gene_results,paste0("/lindstroem/austin_working/Dissertation/AoU/gene_based_results/individual_cancers_combined/",i,"_SAIGE_results.txt"),
#              row.names = FALSE, quote = FALSE, sep = "\t")
#       fwrite(single_results,paste0("/lindstroem/austin_working/Dissertation/AoU/gene_based_results/individual_cancers_combined/",i,"_SAIGE_results.txt.singleAssoc.txt"),
#              row.names = FALSE, quote = FALSE, sep = "\t")
# }

  

# Loop to extract sample sizes and validate results

calculate_lambda_regression <- function(p_values) {
  if (!is.numeric(p_values) || any(p_values <= 0 | p_values > 1)) {
    stop("Input must be a numeric vector of p-values between 0 and 1.")
  }
  
  # Convert p-values to chi-square statistics
  chisq_values <- qchisq(p_values, df = 1, lower.tail = FALSE)
  
  # Sort chi-square values
  chisq_values <- sort(chisq_values)
  
  # Expected chi-square values under the null
  expected_chisq <- qchisq(ppoints(length(chisq_values)), df = 1)
  
  # Perform regression: observed ~ expected, forcing through the origin
  fit <- lm(chisq_values ~ expected_chisq - 1)  # No intercept
  
  # Extract lambda as the slope of the regression line
  lambda <- coef(fit)[1]
  
  return(lambda)
}

validation_data <- data.frame(
  Cancer = character(),  
  N_total = numeric(), 
  N_cases = numeric(), 
  N_controls = numeric(), 
  N_tested_genes = integer(),  
  regression_lambda_synonymous = numeric(),  
  regression_lambda_missense = numeric(),  
  regression_lambda_pLoF = numeric(),  
  stringsAsFactors = FALSE  
)

j = 1
for (i in cancers) {
  gene = fread(paste0("/lindstroem/austin_working/Dissertation/AoU/gene_based_results/individual_cancers_combined/",i,"_SAIGE_results.txt"))
  single = fread(paste0("/lindstroem/austin_working/Dissertation/AoU/gene_based_results/individual_cancers_combined/",i,"_SAIGE_results.txt.singleAssoc.txt"))
  gene_pLoF = gene %>%
    filter(Group == "pLoF" & max_MAF == 0.5)
  
  gene_missense = gene %>%
    filter(Group == "missense" & max_MAF == 0.5)
  
  gene_synonymous = gene %>%
     filter(Group == "synonymous" & max_MAF == 0.01) 
  
  if (j > nrow(validation_data)) {
    validation_data[j, ] <- NA  
  }
  
   validation_data$Cancer[j] <- as.character(cancers[j])
   validation_data$N_total[j] <- median(single$N_case, na.rm = TRUE) + median(single$N_ctrl, na.rm = TRUE)
   validation_data$N_cases[j] <- median(single$N_case, na.rm = TRUE)
   validation_data$N_controls[j] <- median(single$N_ctrl, na.rm = TRUE)
   validation_data$N_tested_genes[j] <- length(unique(gene$Region))
   validation_data$regression_lambda_synonymous[j] <- as.numeric(calculate_lambda_regression(gene_synonymous$Pvalue_Burden))
   validation_data$regression_lambda_missense[j] <- as.numeric(calculate_lambda_regression(gene_missense$Pvalue_Burden))
   validation_data$regression_lambda_pLoF[j] <- as.numeric(calculate_lambda_regression(gene_pLoF$Pvalue_Burden))
   j = j + 1
}

fwrite(validation_data, "/lindstroem/austin_working/Dissertation/AoU/gene_based_results/individual_cancers_combined/UKB_individual_cancer_validation_results.txt",
       row.names = FALSE, quote = FALSE, sep = "\t")
