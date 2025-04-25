#!/bin/bash

# Define paths and samples
VCF="/media/m0hamed/files/files/bioinformatics/wes/ref/trio_annotated_vep.vcf.gz"
PROBAND="proband"
FATHER="father"
MOTHER="mother"
AR_VCF="autosomal_recessive_candidates.vcf"
AR_TSV="autosomal_recessive_candidates.tsv"
CH_TEMP="compound_het_temp.vcf"
CH_TSV="compound_het_candidates.tsv"

# Step 1: Autosomal Recessive
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

# Step 2: Compound Heterozygous
echo "Running compound heterozygous filtering..."
bcftools view \
    -s "$FATHER,$MOTHER,$PROBAND" \
    -i "GT[2]='0/1' && (GT[0]='0/0' || GT[0]='0/1') && (GT[1]='0/0' || GT[1]='0/1') && INFO/CSQ !~ 'IMPACT=LOW'" \
    -o "$CH_TEMP" \
    "$VCF" || {
    echo "Compound heterozygous filtering failed!"
    exit 1
}
bcftools query \
    -f '%CHROM\t%POS\t%REF\t%ALT\t%CSQ\t%ID\tINFO/max_aaf_all\n' \
    "$CH_TEMP" | awk '
    BEGIN {FS="\t"; OFS="\t"}
    {
        split($5, csq, "|"); gene=csq[4];
        if (gene in genes) {
            print genes[gene]; print $0; delete genes[gene];
        } else {
            genes[gene]=$0;
        }
    }' > "$CH_TSV"
echo "Compound heterozygous candidates saved to $CH_TSV"
echo "Top candidates:"
head -n 10 "$CH_TSV" | column -t

echo "Candidate variant detection complete!"
