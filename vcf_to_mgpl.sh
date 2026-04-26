# Unzip the VCF file
gunzip peromyscus_plus_maniculatus_filtered.vcf.gz

# Verify it's now uncompressed
ls -lh peromyscus_plus_maniculatus_filtered.vcf

# Load modules
module load bcftools
module load r/4.3.1

# Create ploidy_inds.txt (if not already exists)
NUM_SAMPLES=$(bcftools query -l peromyscus_plus_maniculatus_filtered.vcf | wc -l)
echo "Number of samples: $NUM_SAMPLES"

# Create ploidy file
for i in $(seq 1 $NUM_SAMPLES); do echo "2"; done > ploidy_inds.txt
echo "Created ploidy_inds.txt with $(wc -l < ploidy_inds.txt) lines"

# Run the R script (the same one that worked for ind995)
Rscript inputdataformat.R peromyscus_plus_maniculatus_filtered.vcf

# Check outputs
echo ""
echo "=== OUTPUT FILES ==="
ls -lh *.mpgl 2>/dev/null
ls -lh pntest_mean_gl.txt 2>/dev/null
ls -lh qk*.txt 2>/dev/null
ls -lh inds_*.txt 2>/dev/null
