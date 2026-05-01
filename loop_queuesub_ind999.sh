#this is a code to submit 3 repeats per entropy run. I'm running entropy for k=2,3 and with 200k, 225k and 250k iterations.
#### noticed that if i loop things every hour instead of every 5 seconds the jobs start faster. For example, if i submit 18 jobs with such high computational demand basically at the same time, the cluster stalls me for 12+ hours, but if i submit them every hour, they start within 30 minutes each
## what i said up here is a LIE!! i submitted the same loop but for ind995 instead of ind999 about 3 days ago and some jobs still haven't started....
#SBATCH --nodes=1
#SBATCH --time=1-00:00:00
#SBATCH --account=def-mcfarlas
#SBATCH --job-name="loop_entropy_ind999"
#SBATCH --ntasks-per-node=1
#SBATCH --mem=4000M
#SBATCH --mail-user=sargon.yousefian@gmail.com
#SBATCH --mail-type=ALL

### ---------- Parameters --------------------------------------------------

MAX_K=3                              # k=2,3 (k=1 causes floating point error)
NUM_REPS=3                           # 3 replicates per k
ITERATIONS=(200000 225000 250000)    # Different iteration lengths
BURNIN=150000                        

# Paths
MPGL_FILE="/home/sydt/scratch/mgpl_file/top_snps_selection/ind999/peromyscus_plus_maniculatus_filtered.mpgl"
STARTING_VALS_DIR="/home/sydt/scratch/mgpl_file/top_snps_selection/ind999"

# Results directory
RESULTS_DIR="/home/sydt/scratch/mgpl_file/top_snps_selection/ind999/entropy_results"
mkdir -p $RESULTS_DIR

echo "========================================="
echo "Starting loop submission script at: $(date)"
echo "Max K value: $MAX_K (running K=2 and K=3 only)"
echo "Number of replicates per K: $NUM_REPS"
echo "Iteration lengths: ${ITERATIONS[@]}"
echo "Burn-in: $BURNIN"
echo "Will submit 1 job per hour"
echo "========================================="

total_jobs=0

# Start from k=2 to skip k=1 (floating point exception)
for k in $(seq 2 $MAX_K)
do
    for iter in "${ITERATIONS[@]}"
    do
        # Convert iterations to readable format
        if [ $iter -eq 200000 ]; then
            ITER_NAME="200k"
        elif [ $iter -eq 225000 ]; then
            ITER_NAME="225k"
        elif [ $iter -eq 250000 ]; then
            ITER_NAME="250k"
        fi
        
        for rep in $(seq 1 $NUM_REPS)
        do 
            JOB_NAME="ind999_k${k}_${ITER_NAME}_rep${rep}"
            echo "Submitting job: $JOB_NAME"
            
            QK_FILE="${STARTING_VALS_DIR}/qk${k}inds.txt"
            
            if [ ! -f "$QK_FILE" ]; then
                echo "WARNING: $QK_FILE not found! Skipping..."
                continue
            fi
            
            sbatch --job-name="$JOB_NAME" \
                   /home/sydt/scratch/mgpl_file/top_snps_selection/ind999/entropy_queuesub_ind999.sh \
                   "$MPGL_FILE" \
                   $k \
                   $rep \
                   "$QK_FILE" \
                   $iter \
                   $BURNIN
            
            ((total_jobs++))
            
            # Wait 1 hour before submitting next job (unless this is the last job)
            if [ $total_jobs -lt 18 ]; then
                echo "Waiting 1 hour before next submission..."
                sleep 3600  # 3600 seconds = 1 hour
            fi
        done
    done
done

echo "Total jobs submitted: $total_jobs"
echo "========================================="
