# Whole Genome Variant Analysis Pipeline

## Introduction
This workflow documents the steps followed to process whole-genome sequencing (WGS) data for a trio (father, mother, and proband). It includes quality control, trimming, alignment, variant calling, annotation, and filtering. Each step includes the rationale behind using specific tools and parameters.

---

## Step 1: Checking Software Versions

### *Commands*
bash
# Check FastQC version
fastqc --version

# Check Trimmomatic version
trimmomatic --version

# Install and check samtools
sudo apt update && sudo apt install samtools -y
samtools --version

# Check BWA version
bwa index --version

# Install Java
apt update && apt install -y openjdk-17-jdk
java -version

# Check GATK version
gatk --version


### *Rationale*
Ensure all tools are installed and compatible with the pipeline.

---

## Step 2: Quality Control (FastQC)

### *Commands*
bash
cd /path/to/raw_data
ls -l *.fq
fastqc *.fq
ls -l *.html *.zip


### *Rationale*
FastQC provides metrics such as base quality, GC content, adapter contamination, and sequence duplication levels.

---

## Step 3: Read Trimming (Trimmomatic)

### *Commands*
bash
trimmomatic PE -phred33 father_R1.fq father_R2.fq -baseout father MINLEN:5
trimmomatic PE -phred33 mother_R1.fq mother_R2.fq -baseout mother MINLEN:5
trimmomatic PE -phred33 proband_R1.fq proband_R2.fq -baseout proband MINLEN:5
ls


### *Rationale*
Removes low-quality bases and adapters. MINLEN:5 ensures very short reads are discarded.

---

## Step 4: Reference Genome Preparation

### *Commands*
bash
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/GRCh38.primary_assembly.genome.fa.gz
gunzip GRCh38.primary_assembly.genome.fa.gz
bwa index GRCh38.primary_assembly.genome.fa
samtools faidx GRCh38.primary_assembly.genome.fa
picard CreateSequenceDictionary R=GRCh38.primary_assembly.genome.fa O=GRCh38.primary_assembly.genome.dict


### *Rationale*
Indexing the genome allows efficient alignment and variant calling.

---

## Step 5: Read Alignment (BWA-MEM)

### *Commands*
bash
bwa mem -t 4 GRCh38.primary_assembly.genome.fa father_R1.fq father_R2.fq > father_aligned.sam
bwa mem -t 4 GRCh38.primary_assembly.genome.fa mother_R1.fq mother_R2.fq > mother_aligned.sam
bwa mem -t 4 GRCh38.primary_assembly.genome.fa proband_R1.fq proband_R2.fq > proband_aligned.sam


### *Rationale*
Aligns reads to the reference genome using BWA-MEM, enabling downstream analysis.

---

## Step 6: Convert SAM to BAM (Samtools View)

### *Commands*
bash
samtools view -S -b father_aligned.sam > father_aligned.bam
samtools view -S -b mother_aligned.sam > mother_aligned.bam
samtools view -S -b proband_aligned.sam > proband_aligned.bam


### *Rationale*
BAM files are more compact and efficient for storage and analysis.

---

## Step 7: Sorting and Indexing BAM Files

### *Commands*
bash
samtools sort father_aligned.bam -o father_sorted.bam
samtools sort mother_aligned.bam -o mother_sorted.bam
samtools sort proband_aligned.bam -o proband_sorted.bam

samtools index father_sorted.bam
samtools index mother_sorted.bam
samtools index proband_sorted.bam

bwa index GRCh38.primary_assembly.genome.fa
ls -lh GRCh38.primary_assembly.genome.fa.*


### *Rationale*
Sorted and indexed BAM files are required for variant calling.

---

## Step 8: Variant Calling (GATK HaplotypeCaller)

### *Commands*
bash
gatk HaplotypeCaller -R GRCh38.primary_assembly.genome.fa -I father_rg.bam -O father.g.vcf.gz -ERC GVCF
gatk HaplotypeCaller -R GRCh38.primary_assembly.genome.fa -I mother_rg.bam -O mother.g.vcf.gz -ERC GVCF
gatk HaplotypeCaller -R GRCh38.primary_assembly.genome.fa -I proband_rg.bam -O proband.g.vcf.gz -ERC GVCF

gatk CombineGVCFs -R GRCh38.primary_assembly.genome.fa --variant father.g.vcf.gz --variant mother.g.vcf.gz --variant proband.g.vcf.gz -O family_combined.g.vcf.gz

gatk GenotypeGVCFs -R GRCh38.primary_assembly.genome.fa -V family_combined.g.vcf.gz -O family_variants.vcf.gz

gatk VariantFiltration -R GRCh38.primary_assembly.genome.fa -V family_variants.vcf.gz -O family_variants_filtered.vcf.gz


### *Rationale*
Generates per-sample variant calls in GVCF format, performs joint genotyping, and applies variant filtration.

---

## Step 9: Variant Annotation and Reporting

### *Commands*
bash
# Install bcftools
sudo apt update
sudo apt install bcftools
bcftools --version

# Download and prepare SnpEff
wget http://downloads.sourceforge.net/project/snpeff/snpEff_latest_core.zip
unzip snpEff_latest_core.zip
cd snpEff
java -jar snpEff.jar

# Extract and annotate proband variants
bcftools view -s proband -Oz -o proband_only.vcf.gz family_variants_filtered.vcf.gz
bcftools index proband_only.vcf.gz

bcftools isec -p isec_output -n=1 proband_only.vcf.gz father.g.vcf.gz mother.g.vcf.gz

java -Xmx4g -jar snpEff.jar -v GRCh38.86 isec_output/0000.vcf > proband_denovo_annotated.vcf

grep "ANN=" proband_denovo_annotated.vcf | head -n 5


### *Filtering Script (Python)*
python
# filter_variants.py

# This script filters variants with HIGH or MODERATE impact
# and specific functional effects (missense, stop gained, etc.)

input_file = "proband_denovo_summary.tsv"
output_file = "denovo_filtered.tsv"

with open(input_file, "r") as infile, open(output_file, "w") as outfile:
    header = infile.readline()
    outfile.write(header)
    
    for line in infile:
        fields = line.strip().split("\t")
        if len(fields) < 8:
            continue
        effect = fields[5]
        impact = fields[6]
        
        if impact in ["HIGH", "MODERATE"]:
            if effect in ["missense_variant", "stop_gained", "frameshift_variant"]:
                outfile.write(line)


### *Run Script*
bash
python filter_variants.py


### *Rationale*
bcftools was used to extract the proband-specific variants and detect de novo ones. SnpEff annotates variants with functional consequences. A custom Python script was used to filter for variants with high or moderate impact, prioritizing biologically significant mutations.

---

## Conclusion

This pipeline successfully processed whole-genome sequencing data for a trio, from raw reads to biologically annotated variants. It integrated industry-standard tools like FastQC, Trimmomatic, BWA, GATK, bcftools, and SnpEff, and provided a systematic approach to variant filtering and annotation.

After annotating the variants, we interpreted the results in a biological context. The focus was placed on non-coding regions, particularly long non-coding RNAs (lncRNAs), due to their emerging role in gene regulation. Most of the identified variants were located in non-coding regions and annotated with a MODIFIER impact level, suggesting potential regulatory functions rather than direct protein-coding alterations.

Among the lncRNAs, *LINC01128* stood out with the highest number of upstream and intronic mutations (42 variants), suggesting a possible role in transcriptional regulation. Recent studies suggest that LINC01128 may be involved in key pathways such as Wnt/β-catenin signaling, which is crucial in bone development and remodeling, indicating a possible connection to osteoporosis.

Other significantly mutated lncRNAs included *LINC00115, **RP5-857K21.4, **RPL7P9, and **NDUFS5P2-RPL7P9*. These lncRNAs showed mutations in intronic and regulatory regions, which may influence alternative splicing or enhancer activity. While RPL7P9 is a pseudogene and NDUFS5P2-RPL7P9 represents a fusion region with unclear function, their mutation patterns suggest potential indirect effects on bone metabolism.

In conclusion, this analysis supports the hypothesis that non-coding variants especially within lncRNAs may contribute to the genetic architecture of complex diseases such as osteoporosis. The results emphasize the importance of expanding genomic studies beyond protein-coding regions to discover novel regulatory elements and mechanisms involved in disease susceptibility and biological function.

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
- [The custom Python script for variant summarization was developed with assistance from ChatGPT]
