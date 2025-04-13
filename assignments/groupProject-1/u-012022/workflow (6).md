
# Whole Genome Variant Analysis Pipeline

## Introduction
This workflow documents the steps followed to process whole-genome sequencing (WGS) data for a trio (father, mother, and proband). It includes quality control, trimming, alignment, variant calling, annotation, and filtering. Each step includes the rationale behind using specific tools and parameters.

## Step 1: Checking Software Versions
Before running the analysis, software versions were checked to ensure compatibility.

### Commands:
```bash
fastqc --version
trimmomatic --version
sudo apt update && sudo apt install samtools -y
samtools --version
bwa index --version
apt update && apt install -y openjdk-17-jdk
java -version
gatk --version
```

### Rationale:
Ensuring that all tools are installed and compatible with the workflow.

## Step 2: Quality Control (FastQC)
FastQC was used to assess the quality of raw sequencing reads.

### Commands:
```bash
cd /path/to/raw_data
ls -l *.fq
fastqc *.fq
ls -l *.html *.zip
```

### Rationale:
FastQC provides metrics such as base quality, GC content, adapter contamination, and sequence duplication levels.

## Step 3: Read Trimming (Trimmomatic)
Trimmomatic was used to remove low-quality reads and adapters.

### Commands:
```bash
trimmomatic PE -phred33 father_R1.fq father_R2.fq -baseout father MINLEN:5
trimmomatic PE -phred33 mother_R1.fq mother_R2.fq -baseout mother MINLEN:5
trimmomatic PE -phred33 proband_R1.fq proband_R2.fq -baseout proband MINLEN:5
ls
```

### Rationale:
Removing low-quality bases improves downstream alignment and variant calling. MINLEN:5 ensures that very short reads are discarded.

## Step 4: Reference Genome Preparation
Downloading and indexing the reference genome (GRCh38).

### Commands:
```bash
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/GRCh38.primary_assembly.genome.fa.gz
gunzip GRCh38.primary_assembly.genome.fa.gz
bwa index GRCh38.primary_assembly.genome.fa
samtools faidx GRCh38.primary_assembly.genome.fa
picard CreateSequenceDictionary R=GRCh38.primary_assembly.genome.fa O=GRCh38.primary_assembly.genome.dict
```

### Rationale:
Indexing the genome allows efficient alignment and variant calling.

## Step 5: Read Alignment (BWA-MEM)
Reads were aligned to the reference genome using BWA-MEM.

### Commands:
```bash
bwa mem -t 4 GRCh38.primary_assembly.genome.fa father_R1.fq father_R2.fq > father_aligned.sam
bwa mem -t 4 GRCh38.primary_assembly.genome.fa mother_R1.fq mother_R2.fq > mother_aligned.sam
bwa mem -t 4 GRCh38.primary_assembly.genome.fa proband_R1.fq proband_R2.fq > proband_aligned.sam
```

### Rationale:
BWA-MEM efficiently maps reads to the reference genome, allowing downstream analysis.

## Step 6: Convert SAM to BAM (Samtools View)
SAM files were converted to BAM format for better storage and processing efficiency.

### Commands:
```bash
samtools view -S -b father_aligned.sam > father_aligned.bam
samtools view -S -b mother_aligned.sam > mother_aligned.bam
samtools view -S -b proband_aligned.sam > proband_aligned.bam
```

### Rationale:
BAM files are more compact and efficient for variant calling.

## Step 7: Sorting and Indexing BAM Files
Sorting ensures alignments are ordered by genomic coordinates.

### Commands:
```bash
samtools sort father_aligned.bam -o father_sorted.bam
samtools sort mother_aligned.bam -o mother_sorted.bam
samtools sort proband_aligned.bam -o proband_sorted.bam
samtools index father_sorted.bam
samtools index mother_sorted.bam
samtools index proband_sorted.bam
bwa index GRCh38.primary_assembly.genome.fa
ls -lh GRCh38.primary_assembly.genome.fa.*
```

### Rationale:
Sorted BAM files are required for downstream variant calling and analysis.

## Step 8: Variant Calling (GATK HaplotypeCaller)
GATK HaplotypeCaller was used to detect variants.

### Commands:
```bash
gatk HaplotypeCaller -R GRCh38.primary_assembly.genome.fa -I father_rg.bam -O father.g.vcf.gz -ERC GVCF
gatk HaplotypeCaller -R GRCh38.primary_assembly.genome.fa -I mother_rg.bam -O mother.g.vcf.gz -ERC GVCF
gatk HaplotypeCaller -R GRCh38.primary_assembly.genome.fa -I proband_rg.bam -O proband.g.vcf.gz -ERC GVCF

gatk CombineGVCFs -R GRCh38.primary_assembly.genome.fa   --variant father.g.vcf.gz   --variant mother.g.vcf.gz   --variant proband.g.vcf.gz   -O family_combined.g.vcf.gz

gatk GenotypeGVCFs -R GRCh38.primary_assembly.genome.fa   -V family_combined.g.vcf.gz   -O family_variants.vcf.gz

gatk VariantFiltration -R GRCh38.primary_assembly.genome.fa   -V family_variants.vcf.gz   -O family_variants_filtered.vcf.gz
```

### Rationale:
Generates per-sample variant calls in GVCF format for joint genotyping, combines GVCF files, performs genotype calling, and applies variant filtration.

## Step 9: Variant Annotation and Reporting

### Commands:
```bash
sudo apt update
sudo apt install bcftools
bcftools --version

wget http://downloads.sourceforge.net/project/snpeff/snpEff_latest_core.zip
unzip snpEff_latest_core.zip
cd snpEff
java -jar snpEff.jar

bcftools view -s proband -Oz -o proband_only.vcf.gz family_variants_filtered.vcf.gz
bcftools index proband_only.vcf.gz

bcftools isec -p isec_output -n=1 proband_only.vcf.gz father.g.vcf.gz mother.g.vcf.gz

java -Xmx4g -jar snpEff.jar -v GRCh38.86 isec_output/0000.vcf > proband_denovo_annotated.vcf

grep "ANN=" proband_denovo_annotated.vcf | head -n 5
```

### Filtering script (filter_variants.py):
```python
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
```

### Run script:
```bash
python filter_variants.py
```

### Rationale:
bcftools was used to extract the proband-specific variants and detect de novo ones. SnpEff annotates variants with functional consequences. A custom Python script was used to filter for variants with high or moderate impact, prioritizing biologically significant mutations.

## Conclusion

This pipeline successfully processed whole-genome sequencing data for a trio, from raw reads to biologically annotated variants. It integrated industry-standard tools like FastQC, Trimmomatic, BWA, GATK, bcftools, and SnpEff, and provided a systematic approach to variant filtering and annotation.

After annotating the variants, we interpreted the results in a biological context. The focus was placed on non-coding regions, particularly long non-coding RNAs (lncRNAs), due to their emerging role in gene regulation. Most of the identified variants were located in non-coding regions and annotated with a MODIFIER  impact level, suggesting potential regulatory functions rather than direct protein-coding alterations.

Among the lncRNAs, LINC01128 stood out with the highest number of upstream and intronic mutations (42 variants), suggesting a possible role in transcriptional regulation. Recent studies suggest that LINC01128 may be involved in key pathways such as Wnt/β-catenin signaling, which is crucial in bone development and remodeling, indicating a possible connection to osteoporosis.

Other significantly mutated lncRNAs included LINC00115, RP5-857K21.4, RPL7P9, and NDUFS5P2-RPL7P9. These lncRNAs showed mutations in intronic and regulatory regions, which may influence alternative splicing or enhancer activity. While RPL7P9 is a pseudogene and NDUFS5P2-RPL7P9 represents a fusion region with unclear function, their mutation patterns suggest potential indirect effects on bone metabolism.

In conclusion, this analysis supports the hypothesis that non-coding variants especially within lncRNAs may contribute to the genetic architecture of complex diseases such as osteoporosis. The results emphasize the importance of expanding genomic studies beyond protein-coding regions to discover novel regulatory elements and mechanisms involved in disease susceptibility and biological function.
