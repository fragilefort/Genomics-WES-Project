#!/bin/bash

# Array of sample names
samples=(father mother proband)
total=${#samples[@]}

# Loop through samples 
for sample in "${samples[@]}"; do
  echo "Processing $sample"
  trimmomatic PE -phred33 \
    ${sample}_R1.fq.gz ${sample}_R2.fq.gz \
    ${sample}_R1_trimmed.fq.gz ${sample}_R1_unpaired.fq.gz \
    ${sample}_R2_trimmed.fq.gz ${sample}_R2_unpaired.fq.gz \
    LEADING:32 TRAILING:32 MINLEN:36 2>&1 | tee ${sample}_trimming.log
  echo "Finished $sample"
done 

echo "All samples trimmed!"
