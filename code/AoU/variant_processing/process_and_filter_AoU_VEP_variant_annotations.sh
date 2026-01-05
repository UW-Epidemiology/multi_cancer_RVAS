#!/bin/bash
#SBATCH --job-name=process_VEP_annotations
#SBATCH --time=20-00:00:00
#SBATCH --cpus-per-task=16
#SBATCH -o logs/%x_%J_out
#SBATCH -e logs/%x_%J_err

Rscript process_and_filter_AoU_VEP_variant_annotations.R