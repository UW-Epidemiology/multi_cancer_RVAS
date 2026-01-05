version 1.0

workflow SAIGE_VR_exome {
  input {
    Array[File] fam_files
    Array[File] bim_files
    Array[File] bed_files
    File sample_file
    Array[Int] chromosome_indices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21] 
    Array[String] instance_mem = ["64GB", "64GB", "32GB", "32GB", "32GB", "32GB", "32GB", "16GB", "32GB", "32GB", "32GB", "32GB", "16GB", "16GB", "16GB", "32GB", "32GB", "16GB", "32GB", "16GB", "16GB", "16GB"]   
    Array[Int] instance_cpu = [16, 16, 8, 8, 8, 8, 8, 4, 8, 8, 8, 8, 4, 4, 4, 8, 8, 4, 8, 4, 4, 4]   
    Array[String] instance_disk = ["local-disk 450 SSD", "local-disk 325 SSD", "local-disk 275 SSD", "local-disk 210 SSD", "local-disk 225 SSD", "local-disk 225 SSD", "local-disk 225 SSD", "local-disk 175 SSD", "local-disk 200 SSD", "local-disk 200 SSD", "local-disk 250 SSD", "local-disk 250 SSD", "local-disk 100 SSD", "local-disk 150 SSD", "local-disk 175 SSD", "local-disk 210 SSD", "local-disk 250 SSD", "local-disk 100 SSD", "local-disk 275 SSD", "local-disk 125 SSD", "local-disk 50 SSD", "local-disk 100 SSD"]  
  }

  scatter (chr_index in chromosome_indices) {
    call plink_thin {
      input:
        fam = fam_files[chr_index],
        bim = bim_files[chr_index],
        bed = bed_files[chr_index],
        samples = sample_file,
        chr = chr_index + 1,
        mem = instance_mem[chr_index],
        cpu = instance_cpu[chr_index],
        disk = instance_disk[chr_index]
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
    String mem
    Int cpu
    String disk
  }

  command <<<
  (
    while true; do 
      echo "---------------------------"
      free -h
      df -h
      echo "---------------------------"
      sleep 600
    done
  ) &
  MONITOR_PID=$! 

    plink2 \
      --bed ~{bed} \
      --bim ~{bim} \
      --fam ~{fam} \
      --keep ~{samples} \
      --no-psam-pheno \
      --thin 0.0005 \
      --make-bed \
      --out c~{chr}_thin

  kill $MONITOR_PID
  >>>

  runtime {
    docker: "austinhamm/saige_v1.4.3_plus_utilities:latest"
    memory: mem
    cpu: cpu
    disks: disk
    bootDiskSizeGb: 10
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
  (
    while true; do 
      echo "---------------------------"
      free -h
      df -h
      echo "---------------------------"
      sleep 600
    done
  ) &
  MONITOR_PID=$! 

    cat ~{write_lines(input_beds)} | sed -e 's/.bed//g' > merge_list.txt
    cat merge_list.txt

    plink \
      --merge-list merge_list.txt \
      --make-bed \
      --new-id-max-allele-len 10000 \
      --out Exome_thin_merge

  kill $MONITOR_PID
  >>>

  runtime {
    docker: "austinhamm/saige_v1.4.3_plus_utilities:latest"
    memory: "32GB"
    cpu: 8
    disks: "local-disk 100 SSD"
    bootDiskSizeGb: 10
  }

  output {
    File merged_bed = "Exome_thin_merge.bed"
    File merged_bim = "Exome_thin_merge.bim"
    File merged_fam = "Exome_thin_merge.fam"
    File merge_log = "Exome_thin_merge.log"
  }
}