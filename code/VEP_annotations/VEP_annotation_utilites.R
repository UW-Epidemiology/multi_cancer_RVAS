library(pacman)
p_load(R.utils, data.table, dplyr)

read_VEP <- function(x, cols = NULL) {
  # Read all lines, remove '##' comment lines, and handle '#' in the header line
  VEP_text <- sub("^#", "", readLines(x)[!startsWith(readLines(x), "##")])
  
  # Load data using fread with specified columns
  fread(text = VEP_text, sep = "\t", na.strings = "-",
        data.table = TRUE, select = cols) 
}


