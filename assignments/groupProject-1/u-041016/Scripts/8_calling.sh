#!/bin/bash

# Absolute path to your working directory
WORKDIR=$(pwd)
REF=$WORKDIR/hg19.fa
SAMPLES=(father mother proband)

# HaplotypeCaller for each sample
for SAMPLE in "${SAMPLES[@]}"; do
  echo "Running HaplotypeCaller for $SAMPLE..."
  docker run -v $WORKDIR:/data -w /data \
    broadinstitute/gatk:latest gatk --java-options "-Xmx6g" HaplotypeCaller \
    -R /data/hg19.fa \
    -I /data/${SAMPLE}_dedup.bam \
    -O /data/${SAMPLE}.g.vcf.gz \
    -ERC GVCF
  echo "HaplotypeCaller complete for $SAMPLE"
done

# Joint Genotyping
echo "Combining gVCFs..."
docker run -v $WORKDIR:/data -w /data \
  broadinstitute/gatk:latest gatk --java-options "-Xmx6g" GenomicsDBImport \
  -R /data/hg19.fa \
  -V /data/father.g.vcf.gz \
  -V /data/mother.g.vcf.gz \
  -V /data/proband.g.vcf.gz \
  --genomicsdb-workspace-path /data/trio_db

echo "Running joint genotyping..."
docker run -v $WORKDIR:/data -w /data \
  broadinstitute/gatk:latest gatk GenotypeGVCFs \
  -R /data/hg19.fa \
  -V gendb:///data/trio_db \
  -O /data/trio.vcf.gz

echo "Variant calling complete!"