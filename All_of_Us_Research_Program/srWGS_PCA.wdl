version 1.0

workflow srWGS_PCA {
  input {
    File fam_file
    File bim_file
    File bed_file
    File sample_file
    File high_ld_regions_hg38
  }

call plink_PCA {
  input:
    fam = fam_file,
    bim = bim_file,
    bed = bed_file,
    samples = sample_file,
    high_ld = high_ld_regions_hg38
}

  output {
    File bed = plink_PCA.LD_bed
    File bim = plink_PCA.LD_bim
    File fam = plink_PCA.LD_fam
    File eigenvalues = plink_PCA.PCA_eigenvalues
    File eigenvectors = plink_PCA.PCA_eigenvectors
  }
}

task plink_PCA {
  input {
    File fam
    File bim
    File bed
    File samples
    File high_ld
  }

  command <<<
    plink2 \
      --bed ~{bed} \
      --bim ~{bim} \
      --fam ~{fam} \
      --keep ~{samples} \
      --chr 1-22 \
      --exclude range ~{high_ld} \
      --geno 0.01 \
      --maf 0.01 \
      --write-snplist \
      --out include

    plink2 \
      --bed ~{bed} \
      --bim ~{bim} \
      --fam ~{fam} \
      --keep ~{samples} \
      --extract include.snplist \
      --indep-pairwise 500 50 0.1 \
      --out LD_pruned

    plink2 \
      --bed ~{bed} \
      --bim ~{bim} \
      --fam ~{fam} \
      --keep ~{samples} \
      --extract LD_pruned.prune.in \
      --make-bed \
      --out LD_pruned

    plink2 \
      --bfile LD_pruned \
      --thin-count 250000 \
      --make-bed \
      --out srWGS_GDA_LD_pruned_and_thinned

    plink2 \
      --bfile srWGS_GDA_LD_pruned_and_thinned \
      --pca 30 approx \
      --out srWGS_PCA
  >>>

  runtime {
    docker: "austinhamm/saige_v1.4.3_plus_utilities:latest"
    memory: "128G"
    cpu: 32
    disks: "local-disk 400 SSD"
    bootDiskSizeGb: "10"
  }

  output {
    File LD_bed = "srWGS_GDA_LD_pruned_and_thinned.bed"
    File LD_bim = "srWGS_GDA_LD_pruned_and_thinned.bim"
    File LD_fam = "srWGS_GDA_LD_pruned_and_thinned.fam"
    File PCA_eigenvalues = "srWGS_PCA.eigenval"
    File PCA_eigenvectors = "srWGS_PCA.eigenvec"
  }
}
