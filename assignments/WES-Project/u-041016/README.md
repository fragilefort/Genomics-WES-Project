# Whole-Exome Sequencing (WES) Analysis Pipeline for Osteopetrosis Variant Detection

This repository contains a pipeline for processing whole-exome sequencing (WES) data from a family trio (proband, father, mother) to identify candidate variants associated with osteopetrosis, a rare genetic bone disorder. The pipeline starts with raw FASTQ files, performs alignment, variant calling, annotation, and filtering, and generates a report of potential causative variants under autosomal recessive and compound heterozygous inheritance models.

## Overview

The pipeline consists of 12 steps, implemented as Bash scripts, to process WES data for a trio where the proband is affected by osteopetrosis, and the consanguineous parents are unaffected. The goal is to identify rare, high-impact variants in known osteopetrosis genes (e.g., `TCIRG1`, `CLCN7`, `OSTM1`, `TNFSF11`, `SNX10`).

## Prerequisites

- **System**: Linux (tested on Ubuntu)

- **Conda**: For managing environments and dependencies

- **Docker**: For running GATK tools

- **Reference Files**:

  - Human genome: `hg19.fa` (e.g., `/media/m0hamed/files/files/bioinformatics/wes/ref/hg19.fa`)
  - VEP cache: GRCh37 (e.g., `/media/m0hamed/files/files/bioinformatics/wes/ref/vep_cache`)

- **Tools** (installed via Conda or Docker):

  - FastQC
  - MultiQC
  - Trimmomatic
  - BWA
  - Samtools
  - Picard
  - GATK (via Docker: `broadinstitute/gatk:latest`)
  - VEP
  - bgzip/tabix
  - vcf2db.py
  - GEMINI
  - bcftools

- **Conda Environments**:

  - `wes_analysis`: For FastQC, MultiQC, Trimmomatic, BWA, Samtools, Picard, VEP, bcftools
  - `vcf2db_env`: For vcf2db.py and GEMINI

  ```bash
  conda create -n wes_analysis -c bioconda fastqc multiqc trimmomatic bwa samtools picard vep bcftools
  conda create -n vcf2db_env python=3.8 gemini vcf2db
  ```

## Pipeline Steps

The pipeline processes FASTQ files through alignment, variant calling, annotation, and variant filtering. Each step is implemented as a Bash script, with quality control performed manually or scripted separately.

### Step 1: Quality Control

- **Tools**: FastQC, MultiQC
- **Description**: Assesses the quality of raw FASTQ files for the trio (`father_R1.fq.gz`, `father_R2.fq.gz`, `mother_R1.fq.gz`, `mother_R2.fq.gz`, `proband_R1.fq.gz`, `proband_R2.fq.gz`) to identify issues like low-quality bases or adapter contamination.
- **Command**:

  ```bash
  mkdir fastqc_results
  fastqc -o fastqc_results *.fq.gz
  multiqc fastqc_results -o multiqc_report
  ```
- **Output**: `fastqc_results/*.html`, `multiqc_report/multiqc_report.html`

### Step 2: Download Data (`1_data.sh`)

- **Description**: Downloads raw FASTQ files for the trio from Zenodo.
- **Code**:

  ```bash
  #!/bin/bash
  wget https://zenodo.org/record/3243160/files/father_R1.fq.gz
  wget https://zenodo.org/record/3243160/files/father_R2.fq.gz
  wget https://zenodo.org/record/3243160/files/mother_R1.fq.gz
  wget https://zenodo.org/record/3243160/files/mother_R2.fq.gz
  wget https://zenodo.org/record/3243160/files/proband_R1.fq.gz
  wget https://zenodo.org/record/3243160/files/proband_R2.fq.gz
  ```
- **Explanation**: Uses `wget` to fetch paired-end FASTQ files for each sample, storing them in the working directory.
- **Output**: `father_R1.fq.gz`, `father_R2.fq.gz`, `mother_R1.fq.gz`, `mother_R2.fq.gz`, `proband_R1.fq.gz`, `proband_R2.fq.gz`
- **Command**:

  ```bash
  ./1_data.sh
  ```

### Step 3: Trim Reads (`2_trim.sh`)

- **Tool**: Trimmomatic
- **Description**: Removes low-quality bases and adapters from FASTQ files to improve alignment quality.
- **Code**:

  ```bash
  #!/bin/bash
  samples=(father mother proband)
  total=${#samples[@]}
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
  ```
- **Explanation**: Loops through samples, running Trimmomatic in paired-end mode to trim reads with quality &lt;32 and length &lt;36, logging results.
- **Output**: `${sample}_R1_trimmed.fq.gz`, `${sample}_R2_trimmed.fq.gz`, `${sample}_R1_unpaired.fq.gz`, `${sample}_R2_unpaired.fq.gz`, `${sample}_trimming.log`
- **Command**:

  ```bash
  ./2_trim.sh
  ```

### Step 4: Align Reads (`3_align.sh`)

- **Tool**: BWA
- **Description**: Aligns trimmed reads to the hg19 reference genome, producing SAM files.
- **Code**:

  ```bash
  #!/bin/bash
  samples=(father mother proband)
  reference=/media/m0hamed/files/files/bioinformatics/wes/ref/hg19.fa
  threads=1
  for sample in "${samples[@]}"; do
    echo "Aligning $sample to hg19..."
    bwa mem -t $threads -R "@RG\tID:${sample}\tSM:${sample}" \
      $reference \
      ${sample}_R1_trimmed.fq.gz ${sample}_R2_trimmed.fq.gz > ${sample}.sam
    echo "Finished aligning $sample"
  done
  echo "Alignment complete for all samples!"
  ```
- **Explanation**: Uses BWA-MEM to align paired-end reads, adding read group metadata (`-R`) and outputting SAM files.
- **Output**: `${sample}.sam`
- **Command**:

  ```bash
  ./3_align.sh
  ```

### Step 5: Filter Reads (`4_filter.sh`)

- **Tool**: Samtools
- **Description**: Filters SAM files to retain paired, mapped reads with high mapping quality (`-q 20`), converting to BAM.
- **Code**:

  ```bash
  #!/bin/bash
  samples=(father mother proband)
  for sample in "${samples[@]}"; do
    echo "Filtering paired reads for $sample..."
    samtools view -b -F 4 -F 8 -q 20 \
      ${sample}.sam > ${sample}_filtered.bam
    echo "Finished filtering $sample"
  done
  echo "Filtering complete for all samples!"
  ```
- **Explanation**: Uses `samtools view` to filter out unmapped (`-F 4`), unpaired (`-F 8`), and low-quality (`-q 20`) reads, producing BAM files.
- **Output**: `${sample}_filtered.bam`
- **Command**:

  ```bash
  ./4_filter.sh
  ```

### Step 6: Sort BAM Files (`5_sort.sh`)

- **Tool**: Samtools
- **Description**: Sorts and indexes filtered BAM files for downstream processing.
- **Code**:

  ```bash
  for sample in father mother proband; do
    samtools sort ${sample}_filtered.bam -o ${sample}_filtered_sorted.bam
    samtools index ${sample}_filtered_sorted.bam
  done
  ```
- **Explanation**: Sorts BAM files by coordinate using `samtools sort` and creates index files with `samtools index`.
- **Output**: `${sample}_filtered_sorted.bam`, `${sample}_filtered_sorted.bam.bai`
- **Command**:

  ```bash
  ./5_sort.sh
  ```

### Step 7: Mark Duplicates (`6_dedup.sh`)

- **Tool**: Picard
- **Description**: Marks and removes duplicate reads from sorted BAM files to reduce PCR artifacts.
- **Code**:

  ```bash
  #!/bin/bash
  samples=(father mother proband)
  for sample in "${samples[@]}"; do
    echo "Processing $sample..."
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
  ```
- **Explanation**: Uses Picard’s `MarkDuplicates` to remove duplicates, generating deduplicated BAMs and metrics files.
- **Output**: `${sample}_dedup.bam`, `${sample}_dedup_metrics.txt`
- **Command**:

  ```bash
  ./6_dedup.sh
  ```

### Step 8: Index BAM Files (`7_index.sh`)

- **Tool**: Samtools
- **Description**: Creates index files for deduplicated BAMs for efficient access.
- **Code**:

  ```bash
  #!/bin/bash
  SAMPLES=(father mother proband)
  for SAMPLE in "${SAMPLES[@]}"; do
    echo "Indexing ${SAMPLE}_dedup.bam..."
    samtools index ${SAMPLE}_dedup.bam
    echo "Indexing complete for $SAMPLE"
  done
  ```
- **Explanation**: Indexes deduplicated BAM files using `samtools index`.
- **Output**: `${sample}_dedup.bam.bai`
- **Command**:

  ```bash
  ./7_index.sh
  ```

### Step 9: Variant Calling (`8_calling.sh`)

- **Tool**: GATK (HaplotypeCaller, GenomicsDBImport, GenotypeGVCFs)
- **Description**: Performs per-sample variant calling, combines gVCFs, and conducts joint genotyping.
- **Code**:

  ```bash
  #!/bin/bash
  WORKDIR=$(pwd)
  REF=$WORKDIR/hg19.fa
  SAMPLES=(father mother proband)
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
  ```
- **Explanation**: Uses GATK’s HaplotypeCaller to generate per-sample gVCFs, combines them with GenomicsDBImport, and performs joint genotyping to produce a trio VCF.
- **Output**: `${sample}.g.vcf.gz`, `trio_db`, `trio.vcf.gz`
- **Command**:

  ```bash
  ./8_calling.sh
  ```

### Step 10: Post-Process Variants (`9_post_proc.sh`)

- **Tool**: GATK (LeftAlignAndTrimVariants)
- **Description**: Splits multiallelic variants and normalizes indels in the joint VCF.
- **Code**:

  ```bash
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
  ```
- **Explanation**: Normalizes indels and splits multiallelic variants using GATK, improving VCF consistency.
- **Output**: `trio_normalized_split.vcf.gz`
- **Command**:

  ```bash
  ./9_post_proc.sh
  ```

### Step 11: Annotate Variants (`10_annotate_VEP.sh`)

- **Tool**: VEP
- **Description**: Annotates variants with functional consequences, ClinVar significance, and population frequencies.
- **Code**:

  ```bash
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
  ```
- **Explanation**: Runs VEP to add annotations (e.g., IMPACT, SYMBOL, CLIN_SIG) to the VCF, then compresses and indexes the output.
- **Output**: `trio_annotated_vep.vcf.gz`, `trio_annotated_vep.vcf.gz.tbi`
- **Command**:

  ```bash
  ./10_annotate_VEP.sh
  ```

### Step 12: Create GEMINI Database (`11_create_db.sh`)

- **Tools**: vcf2db.py, GEMINI
- **Description**: Converts the normalized VCF to a GEMINI SQLite database with pedigree information.
- **Code**:

  ```bash
  #!/bin/bash
  WORKDIR=$(pwd)
  VCF_INPUT="$WORKDIR/trio_normalized_split.vcf.gz"
  PED_FILE="$WORKDIR/trio.ped"
  DB_OUTPUT="$WORKDIR/trio.db"
  echo "Checking VCF file: $VCF_INPUT"
  if [ ! -f "$VCF_INPUT" ]; then
      echo "Error: VCF file $VCF_INPUT not found!"
      exit 1
  fi
  echo "Extracting sample names from VCF..."
  SAMPLES=($(zcat "$VCF_INPUT" | grep "^#CHROM" | cut -f 10-))
  if [ ${#SAMPLES[@]} -ne 3 ]; then
      echo "Error: Expected 3 samples (trio), found ${#SAMPLES[@]}: ${SAMPLES[*]}"
      exit 1
  fi
  PROBAND=${SAMPLES[0]}
  FATHER=${SAMPLES[1]}
  MOTHER=${SAMPLES[2]}
  echo "Detected samples - Proband: $PROBAND, Father: $FATHER, Mother: $MOTHER"
  echo "Creating PED file: $PED_FILE"
  cat > "$PED_FILE" << EOL
  family_id	sample_id	paternal_id	maternal_id	sex	phenotype
  fam1	$PROBAND	$FATHER	$MOTHER	1	2
  fam1	$FATHER	0	0	1	1
  fam1	$MOTHER	0	0	2	1
  EOL
  if [ ! -f "$PED_FILE" ]; then
      echo "Error: Failed to create PED file!"
      exit 1
  fi
  echo "PED file created successfully:"
  cat "$PED_FILE"
  echo "Converting VCF to GEMINI database: $DB_OUTPUT"
  python vcf2db.py \
      "$VCF_INPUT" \
      "$PED_FILE" \
      "$DB_OUTPUT" \
      --expand gt_types --expand gt_depths || {
      echo "vcf2db failed! Check above logs for details."
      exit 1
  }
  echo "Verifying database: $DB_OUTPUT"
  if [ ! -f "$DB_OUTPUT" ]; then
      echo "Error: Database file $DB_OUTPUT not created!"
      exit 1
  fi
  echo "Sample query from database:"
  gemini query -q "SELECT chrom, start, end, ref, alt FROM variants LIMIT 10" "$DB_OUTPUT" || {
      echo "GEMINI query failed! Database might be corrupt."
      exit 1
  }
  echo "Database creation complete! Output: $DB_OUTPUT"
  ```
- **Explanation**: Creates a PED file for the trio and uses `vcf2db.py` to convert the VCF to a GEMINI database, enabling complex variant queries.
- **Output**: `trio.ped`, `trio.db`
- **Command**:

  ```bash
  ./11_create_db.sh
  ```

### Step 13: Find Candidate Variants (`12_find.sh`)

- **Tool**: bcftools
- **Description**: Filters variants for autosomal recessive (proband homozygous, parents heterozygous) and compound heterozygous patterns, focusing on high/moderate impact variants.
- **Code**:

  ```bash
  #!/bin/bash
  VCF="/media/m0hamed/files/files/bioinformatics/wes/ref/trio_annotated_vep.vcf.gz"
  PROBAND="proband"
  FATHER="father"
  MOTHER="mother"
  AR_VCF="autosomal_recessive_candidates.vcf"
  AR_TSV="autosomal_recessive_candidates.tsv"
  CH_TEMP="compound_het_temp.vcf"
  CH_TSV="compound_het_candidates.tsv"
  echo "Running autosomal recessive filtering..."
  bcftools view \
      -s "$FATHER,$MOTHER,$PROBAND" \
      -i "GT[2]='1/1' && GT[0]='0/1' && GT[1]='0/1' && INFO/CSQ !~ 'IMPACT=LOW'" \
      -o "$AR_VCF" \
      "$VCF" || {
      echo "Autosomal recessive filtering failed!"
      exit 1
  }
  bcftools query \
      -f '%CHROM\t%POS\t%REF\t%ALT\t%CSQ\t%ID\tINFO/max_aaf_all\n' \
      "$AR_VCF" > "$AR_TSV"
  echo "Autosomal recessive candidates saved to $AR_TSV"
  echo "Top candidates:"
  head -n 10 "$AR_TSV" | column -t
  echo "Candidate variant detection complete!"
  ```
- **Explanation**: Uses `bcftools` to filter variants based on inheritance patterns, outputting VCF and TSV files. Note: `INFO/max_aaf_all` is incorrect (should use `CSQ`’s `MAX_AF`).
- **Output**: `autosomal_recessive_candidates.vcf`, `autosomal_recessive_candidates.tsv`, `compound_het_temp.vcf`, `compound_het_candidates.tsv`
- **Command**:

  ```bash
  ./12_find.sh
  ```

## File Structure

```
wes_pipeline/
├── 1_data.sh
├── 2_trim.sh
├── 3_align.sh
├── 4_filter.sh
├── 5_sort.sh
├── 6_dedup.sh
├── 7_index.sh
├── 8_calling.sh
├── 9_post_proc.sh
├── 10_annotate_VEP.sh
├── 11_create_db.sh
├── 12_find.sh
├── fastqc_results/
├── multiqc_report/
├── hg19.fa
├── father_R1.fq.gz
├── father_R2.fq.gz
├── mother_R1.fq.gz
├── mother_R2.fq.gz
├── proband_R1.fq.gz
├── proband_R2.fq.gz
├── *.sam
├── *.bam
├── *.vcf.gz
├── trio.ped
├── trio.db
├── autosomal_recessive_candidates.tsv
├── compound_het_candidates.tsv
```

## Usage

1. Clone the repository:

   ```bash
   git clone <repository-url>
   cd wes_pipeline
   ```
2. Set up Conda environments:

   ```bash
   conda env create -n wes_analysis -c bioconda fastqc multiqc trimmomatic bwa samtools picard vep bcftools
   conda env create -n vcf2db_env python=3.8 gemini vcf2db
   ```
3. Ensure `hg19.fa` and VEP cache are in `/media/m0hamed/files/files/bioinformatics/wes/ref/`.
4. Run scripts sequentially:

   ```bash
   bash 1_data.sh
   quality check with fastqc
   bash 2_trim.sh
   # ... continue through 12_find.sh
   ```
5. Perform quality control separately before Step 2.

## Notes

- **Sample Names**: Hardcoded as `father`, `mother`, `proband`. Adjust in scripts if different.
- **Reference Genome**: Uses `hg19.fa`. Update paths in scripts for your setup.
- **VEP Cache**: Assumes GRCh37 cache at `/media/m0hamed/files/files/bioinformatics/wes/ref/vep_cache`.