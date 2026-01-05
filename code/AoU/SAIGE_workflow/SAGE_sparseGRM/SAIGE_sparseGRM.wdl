version 1.0

workflow SAIGE_sparseGRM {
  input {
    File fam_file
    File bim_file
    File bed_file
    File sample_file
  }

  call SAIGE_step0 {
    input:
      fam = fam_file,
      bim = bim_file,
      bed = bed_file,
      samples = sample_file
  }
  
  output {
    File matrix_file = SAIGE_step0.matrix_file
    File matrix_sample_file = SAIGE_step0.matrix_sample_file
  }
}

task SAIGE_step0 {
  input {
    File fam
    File bim
    File bed
    File samples
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

  plink2 \
    --bed ~{bed} \
    --bim ~{bim} \
    --fam ~{fam} \
    --keep ~{samples} \
    --make-bed \
    --out srWGS_EHR

  createSparseGRM.R \
    --bedFile="srWGS_EHR.bed" \
    --bimFile="srWGS_EHR.bim"  \
    --famFile="srWGS_EHR.fam"  \
    --outputPrefix=srWGS_EHR \
    --nThreads=64 \
    --numRandomMarkerforSparseKin=5000 \
    --relatednessCutoff=0.175

  kill $MONITOR_PID
  >>>
  
  runtime {
    docker: "austinhamm/saige_v1.4.3_plus_utilities:latest"
    memory: "240GB"
    cpu: 64
    disks: "local-disk 150 SSD"
    bootDiskSizeGb: 10
  }
  
  output {
    File matrix_file = "srWGS_EHR_relatednessCutoff_0.175_5000_randomMarkersUsed.sparseGRM.mtx"
    File matrix_sample_file = "srWGS_EHR_relatednessCutoff_0.175_5000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt"
  }
}