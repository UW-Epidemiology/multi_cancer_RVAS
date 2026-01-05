# Setup #
source("/projects/lindstroem/UKBB_AoU_WGS_annotations/scripts/VEP_annotation_utilites.R")
library(readr)

library(optparse)

# Set up the command-line argument parser
option_list <- list(
  make_option(c("-c", "--chr"), type = "character", help = "Chromosome value", default = NULL)
)

## Parse the command-line arguments ##
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

## Check if the chromosome value is provided ##
if (is.null(opt$chr)) {
  stop("Chromosome value must be specified using --chr argument")
}

# Define functions #

## Function to return the paths for file parts for each chromosome ##
files_by_pattern <- function(directory, pattern) {
  files <- list.files(path = directory, pattern = pattern, full.names = TRUE, recursive = TRUE)
  return(files)
}

## Function to split variant identifiers into chromosome, position, and alleles ##
split_to_dataframe <- function(x) {
  split_values <- strsplit(x, "-")
  
  df <- do.call(rbind, split_values)
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  
  colnames(df) <- c("chr", "pos", "a1", "a2")
  
  return(df)
}

## Function to combine files by variant position ##
merge_files_by_order <- function(data, file_paths_1, file_paths_2) {
  result <- vector("list", length = nrow(data)) 
  current_file_lines_1 <- NULL
  current_file_lines_2 <- NULL
  current_file_path_1 <- NULL
  current_file_path_2 <- NULL
  column_names <- NULL
  
  for (i in seq_len(nrow(data))) {
    # Select the appropriate file set and path based on source
    if (data$source[i] == "v7") {
      file_path <- file_paths_1[[data$file[i]]]
      current_source_lines <- current_file_lines_1
      current_source_path <- current_file_path_1
      reload_file <- is.null(current_file_path_1) || (file_path != current_file_path_1)
    } else {
      file_path <- file_paths_2[[data$file[i]]]
      current_source_lines <- current_file_lines_2
      current_source_path <- current_file_path_2
      reload_file <- is.null(current_file_path_2) || (file_path != current_file_path_2)
    }
    
    # Load the new file only if needed
    if (reload_file) {
      # Read the gzipped file using readLines and pass to fread via text connection
      file_lines <- readLines(gzfile(file_path))  # Read the compressed file lines
      
      # Filter out metadata lines starting with '##'
      file_lines <- file_lines[!grepl("^##", file_lines)]
      
      # Use text connection to fread for reading tab-delimited data
      file_data <- fread(text = paste(file_lines, collapse = "\n"), sep = "\t", header = TRUE)
      
      # Store the column names from the header and clean up any '#' at the start
      column_names <- sub("^#", "", colnames(file_data))
      
      # Update the current file data for the corresponding source
      if (data$source[i] == "v7") {
        current_file_lines_1 <- file_data
        current_file_path_1 <- file_path
        current_source_lines <- current_file_lines_1
      } else {
        current_file_lines_2 <- file_data
        current_file_path_2 <- file_path
        current_source_lines <- current_file_lines_2
      }
      
      # Debugging: Print when a file is reloaded
      message(paste("Loaded file:", file_path))
    }
    
    # Ensure current_source_lines is not NULL before accessing nrow
    if (!is.null(current_source_lines)) {
      adjusted_line_index <- data$line[i]
      # Ensure adjusted_line_index is within the valid range
      if (adjusted_line_index <= nrow(current_source_lines) && adjusted_line_index > 0) {
        # Select the corresponding row and convert it to a list
        line_fields <- as.list(current_source_lines[adjusted_line_index, ])
        
        # Add the line to result
        result[[i]] <- line_fields
      } else {
        # Handle the case where the line index is out of bounds
        warning(paste("Line index", adjusted_line_index, "is out of bounds for file", file_path))
        result[[i]] <- NULL  # Ensure no null values in result
      }
    }
  }
  
  # Remove NULL values from result
  result <- result[!sapply(result, is.null)]
  
  # Convert the result into a data.frame with appropriate column names
  df_result <- do.call(rbind, result) %>%
    as.data.frame()
  
  # Set appropriate column names
  colnames(df_result) <- column_names
  
  # Return the result as a data.frame
  return(df_result)
}

## Function to split the variants into pieces, merge by line, then write results ## 
apply_function_pieces_and_write <- function(df, num_pieces, file_paths_1, file_paths_2, output_dir, chr) {
  # Calculate the number of rows in each piece
  rows_per_piece <- ceiling(nrow(df) / num_pieces)
  
  # Loop over each piece, slice the data, apply the function, and write to file
  for (i in 1:num_pieces) {
    # Define the range of rows for the current piece
    start_row <- (i - 1) * rows_per_piece + 1
    end_row <- min(i * rows_per_piece, nrow(df))
    
    # Slice the data.frame
    piece_df <- df[start_row:end_row, ]
    
    # Apply the predefined merge_files_by_order function to the piece
    result_df <- merge_files_by_order(piece_df, file_paths_1, file_paths_2)
    
    # Define the output file path for the current piece
    output_file <- file.path(output_dir, paste0("WGS_v8_variants_chr",chr,"_part", i, "_VEP_ANN.tsv.gz"))
    
    # Write the result to a compressed tab-separated file
    fwrite(result_df, file = output_file, sep = "\t",
           quote = FALSE, row.names = FALSE, col.names = TRUE,
                compress = "gzip")
    
    # Print message indicating the file has been written
    cat("Written results for piece", i, "to", output_file, "\n")
  }
}

## Main function ##
merge_annotations = function(chr) {
  
  # Get the list of files for both versions
  chr_v7_files <- files_by_pattern("/projects/lindstroem/UKBB_AoU_WGS_annotations/AoU/block_annotations", 
                                 paste0("WGS_v7.1_variants_chr", chr, "_.*")) %>%
    grep("_warnings\\.txt$", ., invert = TRUE, value = TRUE)

  chr_v8_files <- files_by_pattern("/projects//lindstroem/UKBB_AoU_WGS_annotations/AoU/v8_new_block_annotations", 
                                 paste0("WGS_v8_new_variants_chr", chr, "_.*")) %>%
    grep("_warnings\\.txt$", ., invert = TRUE, value = TRUE)


  # Preallocate lists for identifiers
  chr_v7_identifiers <- vector("list", length(chr_v7_files))
  chr_v8_identifiers <- vector("list", length(chr_v8_files))

  # Process v7 files
  for (i in seq_along(chr_v7_files)) {
    identifiers <- read_VEP(chr_v7_files[i], "Uploaded_variation") %>%
      pull(Uploaded_variation) %>%
      split_to_dataframe() %>%
      mutate(source = "v7", file = i, line = row_number())
  
    chr_v7_identifiers[[i]] <- identifiers
  }

  # Process v8 files
  for (i in seq_along(chr_v8_files)) {
    identifiers <- read_VEP(chr_v8_files[i], "Uploaded_variation") %>%
      pull(Uploaded_variation) %>%
      split_to_dataframe() %>%
      mutate(source = "v8", file = i, line = row_number())
  
    chr_v8_identifiers[[i]] <- identifiers
  }

  # Combine results efficiently
  chr_v7_identifiers <- bind_rows(chr_v7_identifiers)
  chr_v8_identifiers <- bind_rows(chr_v8_identifiers)  
  
  rm(identifiers)
  
  all_chr_identifiers = rbind.data.frame(chr_v7_identifiers,chr_v8_identifiers) %>%
    mutate(pos = as.integer(pos)) %>%
    arrange(pos) %>%
    select(c(source,file,line))
  
  rm(chr_v7_identifiers,chr_v8_identifiers)
  
  apply_function_pieces_and_write(all_chr_identifiers,
                                  num_pieces = 100,
                                  file_paths_1 = chr_v7_files,
                                  file_paths_2 = chr_v8_files,
                                  output_dir = paste0("/projects/lindstroem/UKBB_AoU_WGS_annotations/AoU/v8_final_block_annotations/chr",chr),
                                  chr = chr)
  
}

# Run the function #
## Run the function with the specified chromosome value ##
merge_annotations(chr = opt$chr)
