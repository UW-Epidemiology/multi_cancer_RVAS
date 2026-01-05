#!/bin/bash
#SBATCH --job-name=VEP_ANN_UKBB
#SBATCH --time=20-00:00:00
#SBATCH --cpus-per-task=16
#SBATCH -o logs/%x_%A_%a_out
#SBATCH -e logs/%x_%A_%a_err

# Check if the correct number of arguments are provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <chromosome_value>"
    exit 1
fi

# Get chromosome value from user input
chromosome_value=$1

# Change to the directory where files are located
cd /projects/lindstroem/UKBB_AoU_WGS_annotations/UKBB/WES_coordinates/

# Set input file path
input="WES_variants_chr${chromosome_value}.tsv.gz"

# Set output file path
output="/projects/lindstroem/UKBB_AoU_WGS_annotations/UKBB/WES_annotations/WES_variants_chr${chromosome_value}_VEP_ANN.tsv"

# Print paths for debugging
echo "Input file: $input"
echo "Output file: $output"

# Run VEP
/projects/lindstroem/programs/ensembl-vep/vep \
-i ${input} \
--assembly GRCh38 \
--cache --dir /projects/lindstroem/programs/ensembl-vep/caches/ \
--offline \
--fasta /projects/lindstroem/programs/ensembl-vep/FASTA/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
--pick \
--check_existing \
--show_ref_allele \
--total_length \
--plugin REVEL,file=/projects/lindstroem/programs/ensembl-vep/plugins/REVEL/new_tabbed_revel_grch38.tsv.gz \
--plugin AlphaMissense,file=/projects/lindstroem/programs/ensembl-vep/plugins/AlphaMissense/AlphaMissense_hg38.tsv.gz \
--plugin UTRAnnotator,file=/projects/lindstroem/programs/ensembl-vep/plugins/UTRAnnotator/uORF_5UTR_GRCh38_PUBLIC.txt \
--plugin SpliceVault,file=/projects/lindstroem/programs/ensembl-vep/plugins/SpliceVault/SpliceVault_data_GRCh38.tsv.gz \
--plugin CADD,/projects/lindstroem/programs/CADD-scripts-master/data/prescored/GRCh38_v1.7/no_anno/whole_genome_SNVs.tsv.gz,/projects/lindstroem/programs/CADD-scripts-master/data/prescored/GRCh38_v1.7/no_anno/gnomad.genomes.r4.0.indel.tsv.gz \
--sift b \
--polyphen b \
--variant_class \
--symbol \
--mane \
--numbers \
--protein \
--af \
--af_1kg \
--af_gnomade \
--af_gnomadg \
--max_af \
--dont_skip \
--no_stats \
--fork 16 \
--buffer_size 100000 \
--tab \
-o ${output} 

# Compress the output file
gzip ${output}

