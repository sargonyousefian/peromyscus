#!/bin/bash
#SBATCH --job-name=variant_calling
#SBATCH --account=def-mcfarlas
#SBATCH --time=7-00:00:00                 # 7 days (longer than individual processes)
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --mail-user=sargon.yousefian@gmail.com
#SBATCH --mail-type=END
#SBATCH --output=logs/nextflow_%j.out
#SBATCH --error=logs/nextflow_%j.err

# This is the main Nextflow submission script
# It launches the Nextflow workflow which will submit individual jobs to SLURM

# Create logs and reports directories
mkdir -p logs
mkdir -p /home/sydt/scratch/variant_calling/reports

# Load required modules for Compute Canada
module load nextflow/25.04.6
module load samtools/1.22.1
module load bcftools/1.22

# Set Nextflow environment variables
export NXF_OPTS='-Xms1g -Xmx4g'
export NXF_WORK="/home/sydt/scratch/nextflow_work"

# Create work directory
mkdir -p $NXF_WORK

# Print start information
echo "=========================================="
echo "Variant Calling Pipeline - BCFtools"
echo "=========================================="
echo "Start time: $(date)"
echo ""
echo "Processing ~1500 BAM files from 3 library directories:"
echo "  - /home/sydt/scratch/alignment/library1/aligned"
echo "  - /home/sydt/scratch/alignment/library2/aligned"
echo "  - /home/sydt/scratch/alignment/library3/aligned"
echo ""
echo "Output directory: /home/sydt/scratch/variant_calling"
echo "Work directory: $NXF_WORK"
echo ""
echo "=========================================="
echo ""

# Run the Nextflow pipeline
nextflow run main.nf \
    -profile slurm \
    -work-dir $NXF_WORK \
    -resume

# Check exit status
exit_code=$?

echo ""
echo "=========================================="
echo "Pipeline finished at: $(date)"
echo "Exit code: $exit_code"
if [ $exit_code -eq 0 ]; then
    echo "Status: SUCCESS"
    echo "That's a wrap! Great work on set today guys!"
else
    echo "Status: FAILED"
    echo "Check logs in: logs/"
fi
echo "=========================================="

exit $exit_code
