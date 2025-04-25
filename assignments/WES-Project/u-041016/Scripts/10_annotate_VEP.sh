#!/bin/bash

WORKDIR=$(pwd)
CACHEDIR=/media/m0hamed/files/files/bioinformatics/wes/ref/vep_cache

echo "Running VEP annotation..."
vep --cache \
  --dir_cache $CACHEDIR \
  --assembly GRCh37 \
  --input_file $WORKDIR/trio_normalized_split.vcf.gz \
  --output_file $WORKDIR/trio_annotated_vep.vcf \
  --vcf \
  --fork 1 \
  --everything \
  --verbose \
  --force_overwrite \
  --no_stats \
  --offline \
  || { echo "VEP failed! Check logs above."; exit 1; }

echo "Compressing output..."
bgzip -c $WORKDIR/trio_annotated_vep.vcf > $WORKDIR/trio_annotated_vep.vcf.gz
tabix -p vcf $WORKDIR/trio_annotated_vep.vcf.gz

echo "Annotation and compression complete!"