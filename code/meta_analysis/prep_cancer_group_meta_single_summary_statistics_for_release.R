library(pacman)

p_load(data.table, reactable, R.utils, dplyr)

for (group in c("ectoderm", "mesoderm", "endoderm", "hormone", "infectious", "smoking")) {
  sum_stats <- fread(paste0("/lindstroem/austin_working/Dissertation/Meta/Aim1/cancer_groups/single_meta_", group, ".txt"))
  filtered_sum_stats = sum_stats %>%
    filter(P < 5E-2)
  fwrite(filtered_sum_stats,file = paste0(paste0("/lindstroem/austin_working/Dissertation/Meta/Aim1/cancer_groups_single_for_release/single_meta_", group, "_P_value_under_0.05.txt.gz")),
         compress = "gzip")
}
