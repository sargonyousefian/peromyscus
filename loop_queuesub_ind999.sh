#!/bin/bash

## USAGE: sbatch loop_queuesub_ind999.sh max_k num_reps /path/to/mpgl_and_ldak
## Example: sbatch loop_queuesub_ind999.sh 4 3 /home/sydt/scratch/mgpl_file/ind999
## Tips:
##   - modify k and starting values script/path on command line
##   - do NOT put a / at the end of the path on command line
##   - max_k = highest k value to test (e.g., 4 will run k=1,2,3,4)
##   - num_reps = number of replicates per k (typically 3)

### ---------- Job configuration --------------------------------------------

# Run dependent and permanent parameters
# will be run on complete nodes NOT partial

#SBATCH --nodes=1                       # number of nodes to use            
#SBATCH --time=00-0:15:00 		        # time (DD-HH:MM:SS)
#SBATCH --account=def-mcfarlas          # account name
#SBATCH --job-name="loop_entropy"       # name to display in queue
#SBATCH --ntasks-per-node=1             # tasks per node (one core per node)
#SBATCH --mem=4000M                     # memory per node
#SBATCH --mail-user=sargon.yousefian@gmail.com    # UPDATE THIS: who to email
#SBATCH --mail-type=ALL                 # when to email

### ---------- Main loop ----------------------------------------------------

echo "Starting loop submission script at: $(date)"
echo "Max K value: $1"
echo "Number of replicates per K: $2"
echo "Data directory: $3"
echo "------------------------------------------------"

# Counter for total jobs submitted
total_jobs=0

# Loop through k values from 1 to max_k
for k in $(seq 1 $1)
do
    # Loop through replicates
    for rep in $(seq 1 $2)
    do 
        echo "Submitting job: k=$k, replicate=$rep"
        
        # Submit the job with full paths
        sbatch /home/sydt/scratch/mgpl_file/ind999/entropy_queuesub_ind999.sh \
               $3/*.mpgl \
               $k \
               $rep \
               $3/qk"$k"inds.txt
        
        ((total_jobs++))
    done
done

echo "------------------------------------------------"
echo "Total jobs submitted: $total_jobs"
echo "Loop submission complete at: $(date)"
echo ""
echo "Jobs will output to: /home/sydt/scratch/entropy/ind999"
echo "Check job status with: squeue -u sydt"
