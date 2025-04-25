#!/bin/bash

WORKDIR=$(pwd)

echo "Splitting multiallelics and normalizing indels..."
docker run --memory-swap -1 -v $WORKDIR:/data -w /data \
  broadinstitute/gatk:latest gatk --java-options "-Xmx6g" LeftAlignAndTrimVariants \
  -R /data/hg19.fa \
  -V /data/trio.vcf.gz \
  -O /data/trio_normalized_split.vcf.gz \
  --split-multi-allelics

echo "Normalization and splitting complete!"