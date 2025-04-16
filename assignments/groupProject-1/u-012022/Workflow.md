# Whole Genome Variant Analysis Pipeline
**Author:** Rana Nasser

## Introduction
This workflow documents the steps followed to process whole-genome sequencing (WGS) data for a trio (father, mother, and proband). It includes quality control, trimming, alignment, variant calling, annotation, and filtering. Each step includes the rationale behind using specific tools and parameters.

---

## Step 1: Checking Software Versions
Before running the analysis, software versions were checked to ensure compatibility.

### Commands:
```bash
fastqc --version
java -jar /usr/local/bin/trimmomatic.jar
samtools --version
gatk --version
```

### Rationale:
Ensuring that all tools are installed and compatible with the workflow.

---

## Step 2: Docker Container Setup & Data Binding

```bash
docker run -it -v "$env:USERPROFILE\Desktop\wes:/data" --name wgs_runn_container2 wgs_runn
ls
```

### Rationale:
This command links the local directory to the Docker container for data access.

---

## Step 3: Quality Control (FastQC)

```bash
fastqc /data/*.fq
```

### Rationale:
FastQC provides metrics such as base quality, GC content, adapter contamination, and sequence duplication levels.

---

## Step 4: Read Trimming (Trimmomatic)

```bash
# Father
java -jar /usr/local/bin/trimmomatic.jar PE \
    /data/father_R1.fq /data/father_R2.fq \
    /data/father_R1_paired.fq /data/father_R1_unpaired.fq \
    /data/father_R2_paired.fq /data/father_R2_unpaired.fq \
    ILLUMINACLIP:/usr/local/bin/adapters.fa:2:30:10 \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

# Mother
java -jar /usr/local/bin/trimmomatic.jar PE \
    /data/mother_R1.fq /data/mother_R2.fq \
    /data/mother_R1_paired.fq /data/mother_R1_unpaired.fq \
    /data/mother_R2_paired.fq /data/mother_R2_unpaired.fq \
    ILLUMINACLIP:/usr/local/bin/adapters.fa:2:30:10 \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

# Proband
java -jar /usr/local/bin/trimmomatic.jar PE \
    /data/proband_R1.fq /data/proband_R2.fq \
    /data/proband_R1_paired.fq /data/proband_R1_unpaired.fq \
    /data/proband_R2_paired.fq /data/proband_R2_unpaired.fq \
    ILLUMINACLIP:/usr/local/bin/adapters.fa:2:30:10 \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 \
    -phred33
```

### Rationale:
Removing low-quality bases improves downstream alignment and variant calling.

---

## Step 5: Reference Genome Preparation

```bash
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/GRCh38.primary_assembly.genome.fa.gz
gunzip GRCh38.primary_assembly.genome.fa.gz
bwa index GRCh38.primary_assembly.genome.fa
samtools faidx GRCh38.primary_assembly.genome.fa
picard CreateSequenceDictionary R=GRCh38.primary_assembly.genome.fa O=GRCh38.primary_assembly.genome.dict
```

### Rationale:
Indexing the genome allows efficient alignment and variant calling.

---

## Step 6: Read Alignment (BWA-MEM)

```bash
bwa mem /data/GRCh38.primary_assembly.genome.fa /data/proband_R1_paired.fq /data/proband_R2_paired.fq > /data/proband_aligned.sam
bwa mem /data/GRCh38.primary_assembly.genome.fa /data/father_R1_paired.fq /data/father_R2_paired.fq > /data/father_aligned.sam
bwa mem /data/GRCh38.primary_assembly.genome.fa /data/mother_R1_paired.fq /data/mother_R2_paired.fq > /data/mother_aligned.sam
```

### Rationale:
BWA-MEM efficiently maps reads to the reference genome.

---

## Step 7: Convert SAM to BAM & Sorting

```bash
samtools view -bS /data/proband_aligned.sam > /data/proband_aligned.bam
samtools view -bS /data/father_aligned.sam > /data/father_aligned.bam
samtools view -bS /data/mother_aligned.sam > /data/mother_aligned.bam

samtools sort /data/proband_aligned.bam -o /data/proband_sorted.bam
samtools sort /data/father_aligned.bam -o /data/father_sorted.bam
samtools sort /data/mother_aligned.bam -o /data/mother_sorted.bam

samtools index /data/proband_sorted.bam
samtools index /data/father_sorted.bam
samtools index /data/mother_sorted.bam
```

---

## Step 8: Add Read Groups & Indexing

```bash
gatk AddOrReplaceReadGroups -I /data/proband_sorted.bam -O /data/proband_rg.bam -RGID 1 -RGLB lib1 -RGPL Illumina -RGPU unit1 -RGSM proband
gatk AddOrReplaceReadGroups -I /data/father_sorted.bam -O /data/father_rg.bam -RGID 1 -RGLB lib1 -RGPL Illumina -RGPU unit1 -RGSM father
gatk AddOrReplaceReadGroups -I /data/mother_sorted.bam -O /data/mother_rg.bam -RGID 1 -RGLB lib1 -RGPL Illumina -RGPU unit1 -RGSM mother

samtools index /data/proband_rg.bam
samtools index /data/father_rg.bam
samtools index /data/mother_rg.bam
```

---

## Step 9: Variant Calling (GATK)

```bash
gatk HaplotypeCaller -R GRCh38.primary_assembly.genome.fa -I proband_rg.bam -O proband.g.vcf.gz -ERC GVCF
gatk HaplotypeCaller -R GRCh38.primary_assembly.genome.fa -I father_rg.bam -O father.g.vcf.gz -ERC GVCF
gatk HaplotypeCaller -R GRCh38.primary_assembly.genome.fa -I mother_rg.bam -O mother.g.vcf.gz -ERC GVCF

gatk CombineGVCFs -R GRCh38.primary_assembly.genome.fa \
   --variant father.g.vcf.gz \
   --variant mother.g.vcf.gz \
   --variant proband.g.vcf.gz \
   -O family_combined.g.vcf.gz

gatk GenotypeGVCFs -R GRCh38.primary_assembly.genome.fa -V family_combined.g.vcf.gz -O family_variants.vcf.gz
gatk VariantFiltration -R GRCh38.primary_assembly.genome.fa -V family_variants.vcf.gz -O family_variants_filtered.vcf.gz
```

---

## Step 10: Variant Annotation & Filtering

```bash
bcftools --version
wget http://downloads.sourceforge.net/project/snpeff/snpEff_latest_core.zip
unzip snpEff_latest_core.zip

bcftools view -s proband -Oz -o proband_only.vcf.gz family_variants_filtered.vcf.gz
bcftools index proband_only.vcf.gz
bcftools isec -p isec_output -n=1 proband_only.vcf.gz father.g.vcf.gz mother.g.vcf.gz

java -Xmx4g -jar /data/snpEff/snpEff.jar -v GRCh38.86 isec_output/0000.vcf > proband_denovo_annotated.vcf

grep "ANN=" proband_denovo_annotated.vcf | head -n 5
```

---

## Software Setup Inside Docker (Ubuntu 20.04)

```bash
apt update && apt install -y default-jre bwa samtools bcftools curl unzip

# FastQC
cd /opt
curl -O https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip
unzip fastqc_v0.12.1.zip
chmod +x FastQC/fastqc
ln -s /opt/FastQC/fastqc /usr/local/bin/fastqc

# Trimmomatic
curl -LO http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip
unzip Trimmomatic-0.39.zip
ln -s /opt/Trimmomatic-0.39/trimmomatic-0.39.jar /usr/local/bin/Trimmomatic.jar

# GATK
curl -LO https://github.com/broadinstitute/gatk/releases/download/4.5.0.0/gatk-4.5.0.0.zip
unzip gatk-4.5.0.0.zip
ln -s /opt/gatk-4.5.0.0/gatk /usr/local/bin/gatk

# SnpEff
curl -O https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip
unzip snpEff_latest_core.zip
ln -s /opt/snpEff/snpEff.jar /usr/local/bin/snpEff.jar
```

---

## Conclusion
This pipeline successfully processed whole-genome sequencing data for a trio, from raw reads to biologically annotated variants. It integrated industry-standard tools like FastQC, Trimmomatic, BWA, GATK, bcftools, and SnpEff, and provided a systematic approach to variant filtering and annotation.

After annotating the variants, we interpreted the results in a biological context. The focus was placed on non-coding regions, particularly long non-coding RNAs (lncRNAs), due to their emerging role in gene regulation. Most of the identified variants were located in non-coding regions and annotated with a MODIFIER  impact level, suggesting potential regulatory functions rather than direct protein-coding alterations.

Among the lncRNAs, LINC01128 stood out with the highest number of upstream and intronic mutations (42 variants), suggesting a possible role in transcriptional regulation. Recent studies suggest that LINC01128 may be involved in key pathways such as Wnt/Î²-catenin signaling, which is crucial in bone development and remodeling, indicating a possible connection to osteoporosis.

Other significantly mutated lncRNAs included LINC00115, RP5-857K21.4, RPL7P9, and NDUFS5P2-RPL7P9. These lncRNAs showed mutations in intronic and regulatory regions, which may influence alternative splicing or enhancer activity. While RPL7P9 is a pseudogene and NDUFS5P2-RPL7P9 represents a fusion region with unclear function, their mutation patterns suggest potential indirect effects on bone metabolism.

---

> _For questions or collaborations, feel free to reach out to Rana Nasser._
