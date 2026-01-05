version 1.0

workflow SAIGE_VR_WES {
  input {
    Array[File] fam_files
    Array[File] bim_files
    Array[File] bed_files
    File sample_file
    Array[Int] chromosome_indices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21] 
  }

  scatter (chr_index in chromosome_indices) {
    call plink_thin {
      input:
        fam = fam_files[chr_index],
        bim = bim_files[chr_index],
        bed = bed_files[chr_index],
        samples = sample_file,
        chr = chr_index + 1
    }
  }
  
 call plink_merge {
    input:
      input_beds = plink_thin.output_bed,
      input_bims = plink_thin.output_bim,
      input_fams = plink_thin.output_fam
  }

  output {
    File output_bed = plink_merge.merged_bed
    File output_bim = plink_merge.merged_bim
    File output_fam = plink_merge.merged_fam
  }
}

task plink_thin {
  input {
    File fam
    File bim
    File bed
    File samples
    String chr
  }

  command <<<
    plink2 \
      --bed ~{bed} \
      --bim ~{bim} \
      --fam ~{fam} \
      --keep ~{samples} \
      --no-psam-pheno \
      --thin 0.001 \
      --make-bed \
      --out c~{chr}_thin

  >>>

  runtime {
    docker: "austinhamm/saige_v1.4.3_plus_utilities:latest"
    dx_instance_type: "mem2_ssd2_v2_x8"
    dx_timeout: "10h"
  }

  output {
    File output_bed = "c~{chr}_thin.bed"
    File output_bim = "c~{chr}_thin.bim"
    File output_fam = "c~{chr}_thin.fam"
  }
}

task plink_merge {
  input {
    Array[File] input_beds
    Array[File] input_bims
    Array[File] input_fams
  }

  command <<<
    cat ~{write_lines(input_beds)} | sed -e 's/.bed//g' > merge_list.txt
    cat merge_list.txt

    plink \
      --merge-list merge_list.txt \
      --make-bed \
      --out WES_thin_merge
  >>>

  runtime {
    docker: "austinhamm/saige_v1.4.3_plus_utilities:latest"
    dx_instance_type: "mem2_ssd2_v2_x8"
    dx_timeout: "2h"
  }

  output {
    File merged_bed = "WES_thin_merge.bed"
    File merged_bim = "WES_thin_merge.bim"
    File merged_fam = "WES_thin_merge.fam"
    File merge_log = "WES_thin_merge.log"
  }
}