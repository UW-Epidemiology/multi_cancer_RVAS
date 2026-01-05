#!/bin/bash
#SBATCH --job-name=v8_ann_merge
#SBATCH --time=20-00:00:00
#SBATCH --cpus-per-task=8
#SBATCH -o logs/%x_%A_%a_out
#SBATCH -e logs/%x_%A_%a_err

# Check if the correct number of arguments are provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <chromosome_value>"
    exit 1
fi

# Get chromosome value 
chromosome_value=$1

# Run the R script 
Rscript merge_annotations_for_AoU_v8.R --chr ${1}
