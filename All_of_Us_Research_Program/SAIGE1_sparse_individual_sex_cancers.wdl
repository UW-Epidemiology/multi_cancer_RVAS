version 1.0

workflow SAIGE_sparse_null_models_individual_sex_cancers {
  input {
    File sparseGRM_file
    File sparseGRM_samples_file
    File bed_file
    File bim_file
    File fam_file
    File sample_file
    File Exome_bed_file
    File Exome_bim_file
    File Exome_fam_file
    File pheno_file
    Array[String] phenotypes = ["breast", "prostate", "cervix", "endometrium", "ovary", "testis"]
  }

  call plink_merge {
    input:
      fam = fam_file,
      bim = bim_file,
      bed = bed_file,
      samples = sample_file,
      WES_bed = Exome_bed_file,
      WES_bim = Exome_bim_file,
      WES_fam = Exome_fam_file
  }

  scatter (i in range(length(phenotypes))) {
    call SAIGE_1_sparse_GRM {
      input:
        sparseGRM = sparseGRM_file,
        samples = sparseGRM_samples_file,
        fam = plink_merge.out_fam_file,
        bim = plink_merge.out_bim_file,
        bed = plink_merge.out_bed_file,
        pheno_file = pheno_file,
        pheno = phenotypes[i]
    }
  }

  output {
    Array[File] model_files = flatten(SAIGE_1_sparse_GRM.model_file)  
    Array[File] variance_ratio_files = flatten(SAIGE_1_sparse_GRM.variance_ratio_file)  
  }
}


task plink_merge {
  input {
    File fam
    File bim
    File bed
    File samples
    File WES_fam
    File WES_bim
    File WES_bed
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
      --make-bed \
      --out initial

    plink \
      --bfile "initial" \
      --bmerge ~{WES_bed} ~{WES_bim} ~{WES_fam} \
      --make-bed \
      --out final

  kill $MONITOR_PID
  >>>

  runtime {
    docker: "austinhamm/saige_v1.4.3_plus_utilities:latest"
    memory: "32GB"
    cpu: 8
    disks: "local-disk 100 SSD"
    bootDiskSizeGb: 10
    preemptible: 2
  }

  output {
    File out_fam_file = "final.fam"
    File out_bim_file = "final.bim"
    File out_bed_file = "final.bed" 
  }
}


task SAIGE_1_sparse_GRM {
  input {
    File sparseGRM
    File samples
    File fam
    File bim
    File bed
    File pheno_file
    String pheno
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

    step1_fitNULLGLMM.R --sparseGRMFile=~{sparseGRM} --sparseGRMSampleIDFile=~{samples} --bedFile=~{bed} --bimFile=~{bim} --famFile=~{fam} --phenoFile=~{pheno_file} --phenoCol=~{pheno} --covarColList=age_at_sample,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,PC11,PC12,PC13,PC14,PC15,PC16,PC17,PC18,PC19,PC20,PC21,PC22,PC23,PC24,PC25,PC26,PC27,PC28,PC29,PC30 --sampleIDColinphenoFile=person_id --traitType=binary --outputPrefix=~{pheno}_null_model --useSparseGRMtoFitNULL=TRUE --isCateVarianceRatio=TRUE --useSparseGRMforVarRatio=TRUE --IsOverwriteVarianceRatioFile=TRUE

  kill $MONITOR_PID
  >>>

  runtime {
    docker: "austinhamm/saige_v1.4.3_plus_utilities:latest"
    memory: "140GB"
    cpu: 2
    disks: "local-disk 100 SSD"
    bootDiskSizeGb: 10
    preemptible: 2
  }

  output {
    Array[File] model_file = glob("*.rda")  
    Array[File] variance_ratio_file = glob("*varianceRatio.txt")  
  }
}
