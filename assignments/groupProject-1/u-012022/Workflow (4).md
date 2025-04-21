# Whole Exome Variant Analysis Pipeline

## Introduction
This workflow documents the steps followed to process whole-exome sequencing (WES) data for a trio (father, mother, and proband). It includes quality control, trimming, alignment, variant calling, annotation, and filtering. Each step includes the rationale behind using specific tools and parameters, the expected output, and relevant visualization (if applicable).

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

### Expected Output:
Printed versions of the installed tools.

---

## Step 2: Docker Container Setup & Data Binding

```bash
docker run -it -v "$env:USERPROFILE\Desktop\wes:/data" --name wgs_runn_container2 wgs_runn
ls
```

### Purpose:
Docker Container Setup & Data Binding We run the analysis inside a Docker container to ensure a consistent environment, and we bind the data directory to access files from the host system.

### Rationale:
This command links the local directory to the Docker container for data access.

### Expected Output:
Mounted directory accessible inside the container.

---

## Step 3: Quality Control (FastQC)

```bash
fastqc /data/*.fq
```

### Purpose:
We assess the quality of raw reads to identify any sequencing issues such as low-quality bases or adapter contamination before further processing.

### Rationale:
FastQC provides metrics such as base quality, GC content, adapter contamination, and sequence duplication levels.

### Parameter Explanation:
- This command runs FastQC on all FASTQ files (*.fq) located in the /data/ directory.
- fastqc: A quality control tool used to evaluate the quality of raw sequencing data.
- /data/*.fq: This pattern matches all files in the /data/ directory with the .fq extension, which typically represent raw sequencing reads.

### Expected Output:
HTML and zip reports with quality metrics plots.

### Example Plot: Per Base Sequence Quality ![Proband QC R1](<assignments/groupProject-1/u-012022/plots/fastqc/proband (quality control) R1.png>)

Description:

This plot shows the quality scores across all base positions in the sequencing reads.

- X-axis: Base position in the read (from 1 to ~97 bp).
- Y-axis: Phred quality score:
  - Above 30 (green zone): Very high quality (error probability < 1 in 1000).
  - 20–30 (yellow zone): Moderate quality.
  - Below 20 (red zone): Low quality.

Observations:

- Most of the base positions fall within the green zone, indicating high-quality reads.
- There is a slight drop in quality toward the end of the reads, which is typical in next-generation sequencing data.

Conclusion:

The sequencing data is of good quality and suitable for downstream analysis, such as trimming and alignment.

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

### Purpose:
We trim adapters and low-quality bases from the reads to improve alignment accuracy and overall analysis results.

### Rationale:
Removing low-quality bases improves downstream alignment and variant calling.

### Parameter Explanation:
- PE: Paired-end mode
- ILLUMINACLIP: Removes adapter sequences using a FASTA file (2:30:10 = seed mismatches:palindrome clip threshold:simple clip threshold)
- LEADING:3: Removes low-quality bases from the beginning (quality < 3)
- TRAILING:3: Removes low-quality bases from the end (quality < 3)
- SLIDINGWINDOW:4:15: Scans with a 4-base sliding window, cutting when average quality < 15
- MINLEN:36: Discards reads shorter than 36 bases
- -phred33: Specifies quality encoding (default for modern Illumina machines)

### Expected Output:
Trimmed FASTQ files (paired and unpaired).

---

## Step 5: Reference Genome Preparation

```bash
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/GRCh38.primary_assembly.genome.fa.gz
gunzip GRCh38.primary_assembly.genome.fa.gz
bwa index GRCh38.primary_assembly.genome.fa
samtools faidx GRCh38.primary_assembly.genome.fa
picard CreateSequenceDictionary R=GRCh38.primary_assembly.genome.fa O=GRCh38.primary_assembly.genome.dict
```

### Purpose:
We prepare the reference genome by indexing and creating necessary files so it can be used for efficient and accurate alignment and variant calling.

### Rationale:
Indexing the genome allows efficient alignment and variant calling.

### Parameter Explanation:
- bwa mem: BWA algorithm optimized for longer reads
- reference.fa: Reference genome
- R1/R2_paired.fq: Input reads
- > aligned.sam: Output SAM file

### Expected Output:
Reference genome indexed and ready.

---

## Step 6: Read Alignment (BWA-MEM)

```bash
bwa mem /data/GRCh38.primary_assembly.genome.fa /data/proband_R1_paired.fq /data/proband_R2_paired.fq > /data/proband_aligned.sam
bwa mem /data/GRCh38.primary_assembly.genome.fa /data/father_R1_paired.fq /data/father_R2_paired.fq > /data/father_aligned.sam
bwa mem /data/GRCh38.primary_assembly.genome.fa /data/mother_R1_paired.fq /data/mother_R2_paired.fq > /data/mother_aligned.sam
```

### Purpose:
We align the trimmed reads to the reference genome to determine their genomic positions, which is essential for identifying variants later.

### Rationale:
BWA-MEM efficiently maps reads to the reference genome.

### Parameter Explanation:
- bwa mem: BWA algorithm optimized for longer reads
- reference.fa: Reference genome
- R1/R2_paired.fq: Input reads
- > aligned.sam: Output SAM file

### Expected Output:
SAM files with alignment information.

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

### Purpose:
We convert SAM to BAM because BAM is a compressed binary format that saves space and is required by most downstream tools. We also sort and index the BAM for fast access and visualization.

### Rationale:
Convert SAM to BAM to save space and sort/index for downstream analysis.

### Parameter Explanation:
- view -bS: Converts SAM to BAM
- sort: Sorts BAM file by reference position
- index: Creates BAM index file (.bai)

### Expected Output:
Sorted and indexed BAM files.

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

### Purpose:
We add read group information to the BAM file because GATK and similar tools use it to distinguish between samples and sequencing runs during variant calling.

### Rationale:
Read groups are required for variant calling and downstream processing.

### Parameter Explanation:
- -RGID: Read group ID
- -RGLB: Library
- -RGPL: Platform (e.g., Illumina)
- -RGPU: Platform unit
- -RGSM: Sample name

### Expected Output:
BAM files with read groups and index.

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

### Purpose:
We call genetic variants (SNPs and indels) for each sample, then combine and genotype them together to produce a unified set of variants for the trio.

### Rationale:
GATK tools are used to call, combine, genotype, and filter variants from the trio.

### Parameter Explanation:
- -ERC GVCF: Emits reference confidence model
- --variant: Input GVCF files
- -V: Input VCF
- -O: Output filename

### Expected Output:
Filtered VCF file with called variants.

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

### Purpose:
We annotate the variants to predict their biological effects, and we filter out shared variants to identify potential de novo mutations in the proband.

### Rationale:
- Filter for de novo variants.
- Annotate with functional information using SnpEff.

### Parameter Explanation:
- view -s proband: Selects proband variants
- isec -n=1: Extracts unique variants
- snpEff -v: Runs annotation using specified database (GRCh38.86)

### Expected Output:
Annotated VCF with variant impact prediction.

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
After annotating the variants, we interpreted the results in a biological context. The focus was placed on non-coding regions, particularly long non-coding RNAs (lncRNAs), due to their emerging role in gene regulation. Most of the identified variants were located in non-coding regions and annotated with a MODIFIER impact level, suggesting potential regulatory functions rather than direct protein-coding alterations.

Among the lncRNAs, LINC01128 stood out with the highest number of upstream and intronic mutations (42 variants), suggesting a possible role in transcriptional regulation. Recent studies suggest that LINC01128 may be involved in key pathways such as Wnt/β-catenin signaling, which is crucial in bone development and remodeling, indicating a possible connection to osteoporosis.

Other significantly mutated lncRNAs included LINC00115, RP5-857K21.4, RPL7P9, and NDUFS5P2-RPL7P9. These lncRNAs showed mutations in intronic and regulatory regions, which may influence alternative splicing or enhancer activity. While RPL7P9 is a pseudogene and NDUFS5P2-RPL7P9 represents a fusion region with unclear function, their mutation patterns suggest potential indirect effects on bone metabolism.

---

## References and Tutorials used: 
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)  
- [GATK](https://gatk.broadinstitute.org/)  
- [BWA](http://bio-bwa.sourceforge.net/)  
- [Samtools](http://www.htslib.org/)
- [picard](https://github.com/broadinstitute/picard/releases/latest/download/picard.jar )
- [FastQC & Trimmomatic - Dr. Mohamed Abdelmottaleb's Training]
- [BWA Alignment Tutorial](https://youtu.be/iXFeyexbJ44?si=F9UKO2Y6prLrCeAr)
- [part2](https://youtu.be/P9iF_2nDLGg?si=TF0JQBsmuMR69uQE)  
- [part3](https://youtu.be/XH_GU5bQ7TY?si=4D-YphlqJ_069oTf)
- [GATK](https://youtu.be/iHkiQvxyr5c?si=yyxp7VM_i9QOilj0)
- [Variant Annotation with bcftools & snpEff - Part 1](https://youtu.be/zMNZk14WxXE?si=eLMeYARm-C7jr9Zo)  
- [Variant Annotation with bcftools & snpEff - Part 2](https://youtu.be/-rmreyRAbkE?si=TYOUzU0kEE_vfZ_F)

---

