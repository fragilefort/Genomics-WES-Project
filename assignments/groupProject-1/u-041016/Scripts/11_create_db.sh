#!/bin/bash

# Define working directory
WORKDIR=$(pwd)
VCF_INPUT="$WORKDIR/trio_normalized_split.vcf.gz"
PED_FILE="$WORKDIR/trio.ped"
DB_OUTPUT="$WORKDIR/trio.db"

# Step 2: Verify VCF and extract sample names
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

# Step 3: Create PED file
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

# Step 4: Run vcf2db
echo "Converting VCF to GEMINI database: $DB_OUTPUT"
python vcf2db.py \
    "$VCF_INPUT" \
    "$PED_FILE" \
    "$DB_OUTPUT" \
    --expand gt_types --expand gt_depths || {
    echo "vcf2db failed! Check above logs for details."
    exit 1
}

# Step 5: Verify the database
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