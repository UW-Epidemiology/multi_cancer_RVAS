#!/bin/bash

vep_input_directory="/projects/lindstroem/UKBB_AoU_WGS_annotations/AoU/WGS_v8_new_coordinates/"
vep_output_directory="/projects/lindstroem/UKBB_AoU_WGS_annotations/AoU/v8_new_block_annotations/"

# Allow the user to specify the chromosome value
if [ -z "$1" ]; then
  echo "Usage: $0 <chromosome_value>"
  exit 1
fi
chr_val="$1"

echo "Checking VEP annotations for chromosome ${chr_val}..."

# Step 1: Parse files from the VEP input directory
vep_input_files=()
for filepath in "$vep_input_directory"/*; do
  filename=$(basename "$filepath")
  if [[ "$filename" =~ WGS_v8_new_variants_chr${chr_val}_part([0-9]+)\.tsv\.gz ]]; then
    vep_input_files+=("${BASH_REMATCH[1]}") # Capture the part number
  fi
done

# Step 2: Parse files from the VEP output directory
vep_output_files=()
for filepath in "$vep_output_directory"/*; do
  filename=$(basename "$filepath")
  if [[ "$filename" =~ WGS_v8_new_variants_chr${chr_val}_part([0-9]+)_VEP_ANN\.tsv\.gz ]]; then
    vep_output_files+=("${BASH_REMATCH[1]}") # Capture the part number
  fi
done

# Step 3: Sort the parts numerically
sorted_vep_input_files=($(printf "%s\n" "${vep_input_files[@]}" | sort -n))
sorted_vep_output_files=($(printf "%s\n" "${vep_output_files[@]}" | sort -n))

# Step 4: Determine the range of parts
if [ "${#sorted_vep_input_files[@]}" -gt 0 ]; then
  vep_input_min="${sorted_vep_input_files[0]}"
  vep_input_max="${sorted_vep_input_files[-1]}"
  echo "Range of VEP input files (VEP input directory): ${vep_input_min}-${vep_input_max}"
else
  echo "No VEP input files found in the VEP input directory for chr${chr_val}."
fi

if [ "${#sorted_vep_output_files[@]}" -gt 0 ]; then
  vep_output_min="${sorted_vep_output_files[0]}"
  vep_output_max="${sorted_vep_output_files[-1]}"
  echo "Range of VEP output files (VEP output directory): ${vep_output_min}-${vep_output_max}"
else
  echo "No VEP output files found in the VEP output directory for chr${chr_val}."
fi

# Step 5: Check for missing or extra files
# Check for missing files in the VEP input directory
echo "Checking for missing files in the VEP input directory..."
for part in $(seq "$vep_input_min" "$vep_input_max"); do
  file_path="$vep_input_directory/WGS_v8_new_variants_chr${chr_val}_part${part}.tsv.gz"
  if [ ! -f "$file_path" ]; then
    echo "Missing file: $file_path"
  fi
done

# Check for missing files in the VEP output directory
echo "Checking for missing files in the VEP output directory..."
for part in $(seq "$vep_output_min" "$vep_output_max"); do
  file_path="$vep_output_directory/WGS_v8_new_variants_chr${chr_val}_part${part}_VEP_ANN.tsv.gz"
  if [ ! -f "$file_path" ]; then
    echo "Missing file: $file_path"
  fi
done

# Check for extra files in the VEP input directory
echo "Checking for extra files in the VEP input directory..."
for part in "${sorted_vep_input_files[@]}"; do
  file_path="$vep_input_directory/WGS_v8_new_variants_chr${chr_val}_part${part}.tsv.gz"
  if [[ ! " ${sorted_vep_output_files[@]} " =~ " $part " ]]; then
    echo "Extra file in VEP input directory: $file_path"
  fi
done

# Check for extra files in the VEP output directory
echo "Checking for extra files in the VEP output directory..."
for part in "${sorted_vep_output_files[@]}"; do
  file_path="$vep_output_directory/WGS_v8_new_variants_chr${chr_val}_part${part}_VEP_ANN.tsv.gz"
  if [[ ! " ${sorted_vep_input_files[@]} " =~ " $part " ]]; then
    echo "Extra file in VEP output directory: $file_path"
  fi
done

# Step 6: Report range of file sizes (in MB) for the VEP output directory
if [ "${#sorted_vep_output_files[@]}" -gt 0 ]; then
  file_sizes=()
  for part in "${sorted_vep_output_files[@]}"; do
    file_path="$vep_output_directory/WGS_v8_new_variants_chr${chr_val}_part${part}_VEP_ANN.tsv.gz"
    file_size=$(stat -c%s "$file_path")
    file_sizes+=("$file_size")
  done
  
  # Convert file sizes to MB
  file_sizes_mb=()
  for size in "${file_sizes[@]}"; do
    size_in_mb=$(echo "$size / 1048576" | bc -l) # Convert bytes to MB
    file_sizes_mb+=($(printf "%.2f" "$size_in_mb"))
  done

  # Sort file sizes in MB
  sorted_file_sizes_mb=($(printf "%s\n" "${file_sizes_mb[@]}" | sort -n))

  # Calculate range of file sizes
  min_size="${sorted_file_sizes_mb[0]}"
  max_size="${sorted_file_sizes_mb[-1]}"

  # Report the range for file sizes
  echo "Range of VEP output file sizes: ${min_size} MB - ${max_size} MB"
  
  # Separate report for the last part (final file)
  final_part="${sorted_vep_output_files[-1]}"
  final_file_path="$vep_output_directory/WGS_v8_new_variants_chr${chr_val}_part${final_part}_VEP_ANN.tsv.gz"
  final_file_size=$(stat -c%s "$final_file_path")
  final_file_size_mb=$(echo "$final_file_size / 1048576" | bc -l)
  final_file_size_mb_rounded=$(printf "%.2f\n" "$final_file_size_mb")
  echo "Final VEP output file size (part ${final_part}): ${final_file_size_mb_rounded} MB"


  # Calculate the median file size
  sorted_file_sizes_mb=($(printf "%s\n" "${file_sizes_mb[@]}" | sort -n))
  num_files=${#sorted_file_sizes_mb[@]}

  if [ $((num_files % 2)) -eq 1 ]; then
    # Odd number of elements, median is the middle element
    median_size=${sorted_file_sizes_mb[$((num_files / 2))]}
  else
    # Even number of elements, median is the average of the two middle elements
    mid1=$((num_files / 2 - 1))
    mid2=$((num_files / 2))
    median_size=$(echo "(${sorted_file_sizes_mb[$mid1]} + ${sorted_file_sizes_mb[$mid2]}) / 2" | bc -l)
  fi
  median_size_rounded=$(printf "%.2f" "$median_size")
  echo "Median file size: ${median_size_rounded} MB"

  # Calculate 10% deviation from the median
  deviation_threshold=$(echo "$median_size * 0.2" | bc -l)
  lower_threshold=$(echo "$median_size - $deviation_threshold" | bc -l)
  upper_threshold=$(echo "$median_size + $deviation_threshold" | bc -l)

  # List files whose sizes deviate by more than 20% from the median file size
  echo "Files with sizes deviating by more than 20% from the median file size:"
  for i in "${!file_sizes_mb[@]}"; do
    file_size_mb="${file_sizes_mb[$i]}"
    part="${sorted_vep_output_files[$i]}"
    if (( $(echo "$file_size_mb < $lower_threshold" | bc -l) )) || (( $(echo "$file_size_mb > $upper_threshold" | bc -l) )); then
      echo "Part ${part}: ${file_size_mb} MB (deviates from median file size)"
    fi
  done
else
  echo "No VEP output files found to calculate file size range."
fi


