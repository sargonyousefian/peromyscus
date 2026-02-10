#!/bin/bash

# Process first Peromyscus VCF file to MGPL format
# File: peromyscus_mac3_Q30_miss0.25_dp3_ind995_maf001.recode.vcf
# Using K=1,2,3,4 version

set -e  # Exit on error
module load r/4.3.1

VCF_FILE="peromyscus_mac3_Q30_miss0.25_dp3_ind999_maf001.recode.vcf"
PLOIDY=2

echo "=== Processing VCF to MGPL Conversion ==="
echo "File: $VCF_FILE"
echo "Ploidy: $PLOIDY"
echo "K values: 1, 2, 3, 4"
echo ""

# Check if VCF file exists
if [ ! -f "$VCF_FILE" ]; then
    echo "ERROR: VCF file not found: $VCF_FILE"
    echo "Please make sure the file is in the current directory"
    exit 1
fi

# Step 1: Extract sample names and create ploidy file
echo "Step 1: Creating ploidy_inds.txt file..."
# Try bcftools first, fall back to grep if not available
if command -v bcftools &> /dev/null; then
    bcftools query -l "$VCF_FILE" > sample_names.tmp
else
    # Extract sample names from #CHROM line
    grep "^#CHROM" "$VCF_FILE" | cut -f10- | tr '\t' '\n' > sample_names.tmp
fi

NUM_SAMPLES=$(wc -l < sample_names.tmp)
echo "Found $NUM_SAMPLES individuals in VCF file"

# Create ploidy file (one ploidy value per individual)
for i in $(seq 1 $NUM_SAMPLES); do
    echo "$PLOIDY"
done > ploidy_inds.txt

echo "Created ploidy_inds.txt with ploidy=$PLOIDY for all $NUM_SAMPLES individuals"
rm sample_names.tmp

# Step 2: Run R script to convert VCF to MGPL
echo ""
echo "Step 2: Running R script to convert VCF to MGPL format..."
echo "This may take a few minutes depending on file size..."
echo ""

Rscript inputdataformat.R "$VCF_FILE"

echo ""
echo "=== Conversion Complete! ==="
echo ""
echo "Output files created:"
ls -lh *.mpgl 2>/dev/null || echo "  - MGPL file (check for .mpgl extension)"
ls -lh pntest_mean_gl.txt 2>/dev/null || true
ls -lh qk*.txt 2>/dev/null || true
ls -lh inds_*.txt 2>/dev/null || true
echo ""
echo "Main output: peromyscus_mac3_Q30_miss0.25_dp3_ind995_maf001.recode.mpgl"
echo "K values generated: qk1inds.txt, qk2inds.txt, qk3inds.txt, qk4inds.txt"
