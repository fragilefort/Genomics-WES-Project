#!/bin/bash

# Activate Conda environment (run manually if not in script)
# conda activate wes_analysis

# Array of sample names
samples=(father mother proband)

# Path to hg19 reference (adjust if in a different directory)
reference=/media/m0hamed/files/files/bioinformatics/wes/ref/hg19.fa

# Number of threads (adjust based on your CPU)
threads=1

# Loop through samples
for sample in "${samples[@]}"; do
  echo "Aligning $sample to hg19..."
  bwa mem -t $threads -R "@RG\tID:${sample}\tSM:${sample}" \
    $reference \
    ${sample}_R1_trimmed.fq.gz ${sample}_R2_trimmed.fq.gz > ${sample}.sam
  echo "Finished aligning $sample"
done

echo "Alignment complete for all samples!"
