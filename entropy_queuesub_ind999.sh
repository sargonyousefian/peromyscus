#!/bin/bash

#SBATCH --nodes=1
#SBATCH --time=7-00:00:00
#SBATCH --account=def-mcfarlas
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --mail-user=sargon.yousefian@gmail.com
#SBATCH --mail-type=ALL

# Load required modules
module load StdEnv/2020 gcc/9.3.0 gsl/2.6 hdf5/1.10.6

### ---------- Useful job information -------------------------------------

echo "Current working directory: $PWD"
echo "starting run at: $(date)"
echo "------------------------------------------------"
echo "job is running on node: $HOSTNAME"
echo "------------------------------------------------"
echo "Job identifier is $SLURM_JOB_ID"
echo "Job name is $SLURM_JOB_NAME"
echo "Using $SLURM_CPUS_PER_TASK CPU cores"
echo "Allocated memory: $SLURM_MEM_PER_NODE MB"

### ---------- Main ---------------------------------------------------------

# Parse arguments
IN_FILE="$1"
k="$2"
rep="$3"
LDAK_FILE="$4"
ITERATIONS="$5"
BURNIN="$6"

# Debug output
echo "Debug: IN_FILE = '$IN_FILE'"
echo "Debug: k = '$k'"
echo "Debug: rep = '$rep'"
echo "Debug: LDAK_FILE = '$LDAK_FILE'"
echo "Debug: ITERATIONS = '$ITERATIONS'"
echo "Debug: BURNIN = '$BURNIN'"

# Set working directory for results
RESULTS_DIR="/home/sydt/scratch/mgpl_file/top_snps_selection/ind999/entropy_results"
cd $RESULTS_DIR

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
ENTROPY="/home/sydt/projects/def-mcfarlas/software/entropy"

# Extract prefix
PREFIX=$(basename "$IN_FILE" .mpgl)
echo "Debug: PREFIX = '$PREFIX'"

# Convert iterations to readable format
if [ "$ITERATIONS" -eq 200000 ]; then 
    ITER_READABLE="200k"
elif [ "$ITERATIONS" -eq 225000 ]; then 
    ITER_READABLE="225k"
elif [ "$ITERATIONS" -eq 250000 ]; then 
    ITER_READABLE="250k"
else
    ITER_READABLE="${ITERATIONS}"
fi

# Output files
OUTPUT_FILE="${PREFIX}_k${k}_${ITER_READABLE}_rep${rep}.hdf5"
DIC_FILE="${PREFIX}_k${k}_${ITER_READABLE}_rep${rep}_DIC.txt"

echo "Starting entropy with:"
echo "  k = $k, replicate = $rep"
echo "  iterations = $ITERATIONS ($ITER_READABLE)"
echo "  burnin = $BURNIN"
echo "  ldak file = $LDAK_FILE"
echo "Output HDF5: $OUTPUT_FILE"
echo "Output DIC: $DIC_FILE"
echo "------------------------------------------------"

# Run entropy
$ENTROPY -i "$IN_FILE" \
         -l "$ITERATIONS" \
         -b "$BURNIN" \
         -t 25 \
         -k "$k" \
         -o "$OUTPUT_FILE" \
         -n 2 \
         -m 1 \
         -w 1 \
         -q "$LDAK_FILE" \
         -Q 1 \
         -D 1 \
         -r $RANDOM > "$DIC_FILE" 2>&1

# Extract DIC values
grep -E "DIC|deviance|pD" "$DIC_FILE" > "${DIC_FILE}.clean" 2>/dev/null

echo "entropy run complete at: $(date)"
echo "Results in: $RESULTS_DIR"
ls -lh "$OUTPUT_FILE" 2>/dev/null || echo "WARNING: Output file not created!"
echo "------------------------------------------------"
