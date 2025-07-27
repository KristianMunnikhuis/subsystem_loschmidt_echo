#!/bin/bash -l

#$ -P fheating
#$ -m ea
#$ -N Equilibrium_Correlations
#$ -j y

echo "================"
echo "Start date: $(date)"
echo "Job name : $JOB_NAME"
echo "Job ID : $JOB_ID  $SGE_TASK_ID"
echo "=========================================================="

n_values=(8 9 10 11 12)  
for n in "${n_values[@]}"; do
    echo "Running with n = $n"
    python equilibrium_correlations.py "$n"
done

echo "Completed!"

