#!/bin/bash

## Bash script for submitting a job to the Compute Canada queue
## Usage: sbatch entropy_queuesub_ind999.sh mpgl_file k_number rep_number ldak_file

### ---------- Job configuration --------------------------------------------

#SBATCH --nodes=1
#SBATCH --time=7-00:00:00
#SBATCH --account=def-mcfarlas
#SBATCH --job-name="entropy_ind999"
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --mail-user=sargon.yousefian@gmail.com
#SBATCH --mail-type=ALL

# Load required modules
module load nixpkgs/16.09  intel/2017.1 gsl/2.3 hdf5/1.8.18

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

# Set working directory to output location
WORKDIR=/home/sydt/scratch/entropy/ind999
cd $WORKDIR

# Enable OpenMP threading for parallel operations
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Path to entropy program
ENTROPY="/home/sydt/projects/def-mcfarlas/software/entropy"

# Path to mpgl file
IN_FILE=$1
echo "Path to MPGL file: $IN_FILE"
PREFIX=$(basename $IN_FILE .recode.mpgl)
echo "MPGL file without path: $PREFIX"

# Path to ldak file (qk file)
LDAK_FILE=$4
echo "Path to LDAK file: $LDAK_FILE"
SUFFIX=$(basename $LDAK_FILE .txt)
echo "LDAK file without path: $SUFFIX"

echo "starting entropy"

# Entropy call
k=$2
rep=$3

# Print command
entropyrun="$ENTROPY -i $IN_FILE -l 150000 -b 125000 -t 25 -k $k -o ${PREFIX}_k${k}_150k_rep${rep}_${SUFFIX}.hdf5 -n 2 -m 1 -w 0 -q $LDAK_FILE -Q 0 -r $RANDOM"
echo $entropyrun

# Execute command
$ENTROPY -i $IN_FILE -l 150000 -b 125000 -t 25 -k $k -o ${PREFIX}_k${k}_150k_rep${rep}_${SUFFIX}.hdf5 -n 2 -m 1 -w 0 -q $LDAK_FILE -Q 0 -r $RANDOM

echo "entropy run done. Results in ${PREFIX}_k${k}_150k_rep${rep}_${SUFFIX}.hdf5"
echo "Output saved to: $WORKDIR"
echo "Finished at: $(date)"
