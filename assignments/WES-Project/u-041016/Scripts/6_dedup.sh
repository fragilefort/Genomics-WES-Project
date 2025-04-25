#!/bin/bash

# Activate Conda environment (run manually if not in script)
# conda activate wes_analysis

# Array of sample names
samples=(father mother proband)

# Loop through samples
for sample in "${samples[@]}"; do
  echo "Processing $sample..."

  # Mark duplicates
  echo "  Marking duplicates for $sample..."
  picard MarkDuplicates \
    INPUT=${sample}_filtered_sorted.bam \
    OUTPUT=${sample}_dedup.bam \
    METRICS_FILE=${sample}_dedup_metrics.txt \
    REMOVE_DUPLICATES=true \
    CREATE_INDEX=true

  echo "Finished processing $sample"
done

echo "Duplicate marking complete for all samples!"
