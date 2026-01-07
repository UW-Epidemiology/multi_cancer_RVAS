version 1.1

workflow SAIGE1_ancestry_stratified_breast_prostate {
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
    Array[String] phenotypes
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

  scatter (pheno in phenotypes) {
    String mem =
      if (pheno.endsWith("EUR")) {
        "64GB"
      } else if (pheno.endsWith("AFR") || pheno.endsWith("AMR")) {
        "32GB"
      } else {
        "16GB"
      }

    call SAIGE_1_sparse_GRM {
      input:
        sparseGRM = sparseGRM_file,
        samples = sparseGRM_samples_file,
        fam = plink_merge.out_fam_file,
        bim = plink_merge.out_bim_file,
        bed = plink_merge.out_bed_file,
        pheno_file = pheno_file,
        pheno = pheno,
        memory = mem
    }
  }

  output {
    Array[File] model_files = flatten(SAIGE_1_sparse_GRM.model_file)
    Array[File] variance_ratio_files = flatten(SAIGE_1_sparse_GRM.variance_ratio_file)
  }
}
