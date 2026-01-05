version 1.0

workflow SAIGE_GENE_sparse_individual_cancers {
  input {
    Array[File] bed_files
    Array[File] bim_files
    Array[File] fam_files
    Array[File] group_files
    Array[File] variant_lists
    Array[File] model_files
    Array[File] VR_files
    File sparseGRM_file
    File sparseGRM_samples_file
    Array[String] phenotypes = ["anus", "bladder", "bone", "breast", "brain", "cervix", "colorectal", "endometrium", "esophagus", "HL", "kidney", "leukemias", "liver", "lung", "melanoma", "myeloma", "NHL", "neck", "oral", "ovary", "pancreas", "prostate", "stomach", "testis", "thyroid"]
    Array[Int] chromosome_indices = [0, 1, 2, 3, 5, 6, 8, 10, 11, 12, 14, 15, 16, 18, 19, 21] 
    Array[String] instance_mem = ["64GB", "64GB", "32GB", "32GB", "32GB", "32GB", "32GB", "16GB", "32GB", "32GB", "32GB", "32GB", "16GB", "16GB", "16GB", "32GB", "32GB", "16GB", "32GB", "16GB", "16GB", "16GB"]   
    Array[Int] instance_cpu = [16, 16, 8, 8, 8, 8, 8, 4, 8, 8, 8, 8, 4, 4, 4, 8, 8, 4, 8, 4, 4, 4]   
    Array[String] instance_disk1 = ["local-disk 500 SSD", "local-disk 450 SSD", "local-disk 400 SSD", "local-disk 350 SSD", "local-disk 450 SSD", "local-disk 450 SSD", "local-disk 450 SSD", "local-disk 250 SSD", "local-disk 350 SSD", "local-disk 350 SSD", "local-disk 450 SSD", "local-disk 450 SSD", "local-disk 150 SSD", "local-disk 250 SSD", "local-disk 300 SSD", "local-disk 350 SSD", "local-disk 450 SSD", "local-disk 100 SSD", "local-disk 450 SSD", "local-disk 200 SSD", "local-disk 100 SSD", "local-disk 200 SSD"]
    Array[String] instance_disk2 = ["local-disk 225 SSD", "local-disk 175 SSD", "local-disk 135 SSD", "local-disk 105 SSD", "local-disk 113 SSD", "local-disk 113 SSD", "local-disk 113 SSD", "local-disk 100 SSD", "local-disk 100 SSD", "local-disk 100 SSD", "local-disk 125 SSD", "local-disk 125 SSD", "local-disk 100 SSD", "local-disk 100 SSD", "local-disk 100 SSD", "local-disk 105 SSD", "local-disk 125 SSD", "local-disk 100 SSD", "local-disk 135 SSD", "local-disk 100 SSD", "local-disk 100 SSD", "local-disk 100 SSD"]  
  }
  

  scatter (chr_index in chromosome_indices) {
    call plink_extract {
      input:
        fam = fam_files[chr_index],
        bim = bim_files[chr_index],
        bed = bed_files[chr_index],
        variants = variant_lists[chr_index],
        chr = chr_index + 1,
        mem = instance_mem[chr_index],
        cpu = instance_cpu[chr_index],
        disk = instance_disk1[chr_index]
    }

    scatter (pheno_index in range(length(phenotypes))) {
      call SAIGE_2 {
        input:
          bed = plink_extract.out_bed,
          bim = plink_extract.out_bim,
          fam = plink_extract.out_fam,
          phenotype = phenotypes[pheno_index],
          sparseGRM = sparseGRM_file,
          samples = sparseGRM_samples_file,
          model = model_files[pheno_index],
          VR = VR_files[pheno_index],
          group = group_files[chr_index], 
          chr = chr_index + 1, 
          disk = instance_disk2[chr_index]

      }
    }
  }

  output {
    Array[Array[File]] output_txt = SAIGE_2.results_txt
    Array[Array[File]] output_txt_markers = SAIGE_2.results_txt_markers
    Array[Array[File]] output_txt_single = SAIGE_2.results_single_txt
  }
}

task plink_extract {
  input {
    File fam
    File bim
    File bed
    File variants
    Int chr
    String mem
    Int cpu
    String disk
 
  }

  command <<<
    plink2 \
      --bed ~{bed} \
      --bim ~{bim} \
      --fam ~{fam} \
      --no-psam-pheno \
      --extract ~{variants} \
      --set-all-var-ids chr@:#:\$r:\$a \
      --new-id-max-allele-len 10000 \
      --make-bed \
      --out c~{chr}

  >>>

  runtime {
    docker: "austinhamm/saige_v1.4.3_plus_utilities:latest"
    memory: mem
    cpu: cpu
    disks: disk
    bootDiskSizeGb: 10
    preemptible: 2
  }

  output {
    File out_bed = "c~{chr}.bed"
    File out_bim = "c~{chr}.bim"
    File out_fam = "c~{chr}.fam"
  }
}

task SAIGE_2 {
  input {
    File fam
    File bim
    File bed
    String phenotype
    File sparseGRM
    File samples
    File model
    File VR
    File group
    Int chr
    String disk
  }

  command <<<
  (
    while true; do 
      echo "---------------------------"
      free -h
      df -h
      echo "---------------------------"
      sleep 1200
    done
  ) &
    
  MONITOR_PID=$! 

step2_SPAtests.R --bedFile=~{bed} --bimFile=~{bim} --famFile=~{fam} --AlleleOrder=alt-first --SAIGEOutputFile="chr~{chr}_~{phenotype}_SAIGE_results.txt" --chrom=~{chr} --minMAF=0 --minMAC=0.5 --GMMATmodelFile=~{model} --varianceRatioFile=~{VR} --sparseGRMFile=~{sparseGRM} --sparseGRMSampleIDFile=~{samples} --is_Firth_beta=TRUE --groupFile=~{group} --annotation_in_groupTest="pLoF,missense,synonymous,pLoF:missense,missense:synonymous,pLoF:missense:synonymous:other" --maxMAF_in_groupTest=0.01,0.5 --is_output_markerList_in_groupTest=TRUE --is_single_in_groupTest=TRUE --is_output_moreDetails=TRUE --r.corr=1 --LOCO=FALSE --is_fastTest=TRUE

  kill $MONITOR_PID
  >>>

  runtime {
    docker: "wzhou88/saige:1.3.6"
    memory: "64GB"
    cpu: 2
    disks: disk
    bootDiskSizeGb: 10
    preemptible: 2
  }

  output {
    File results_txt = "chr~{chr}_~{phenotype}_SAIGE_results.txt"
    File results_txt_markers = "chr~{chr}_~{phenotype}_SAIGE_results.txt.markerList.txt"
    File results_single_txt = "chr~{chr}_~{phenotype}_SAIGE_results.txt.singleAssoc.txt"
  }
}