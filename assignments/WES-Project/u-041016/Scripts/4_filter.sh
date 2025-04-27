#!/bin/bash

# Activate Conda environment (run manually if not in script)
# conda activate wes_analysis

# Array of sample names
samples=(father mother proband)

# Loop through samples
for sample in "${samples[@]}"; do
  echo "Filtering paired reads for $sample..."
  
  # Filter SAM to BAM with paired, mapped reads
  samtools view -b -F 4 -F 8 -q 20 \
    ${sample}.sam > ${sample}_filtered.bam
  
  echo "Finished filtering $sample"
done

echo "Filtering complete for all samples!"
