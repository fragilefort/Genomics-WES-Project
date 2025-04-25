#!/bin/bash

# Array of sample names
SAMPLES=(father mother proband)

# Index each BAM file
for SAMPLE in "${SAMPLES[@]}"; do
  echo "Indexing ${SAMPLE}_dedup.bam..."
  samtools index ${SAMPLE}_dedup.bam
  echo "Indexing complete for $SAMPLE"
done