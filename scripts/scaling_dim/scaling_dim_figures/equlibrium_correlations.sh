#!/bin/bash -l

#$ -P fheating
#$ -m ea
#$ -N Equilbrium Correlations
#$ -j y


echo "================"
echo "Start date: $(date)"
echo "Job name : $JOB_NAME"
echo "Job ID : $JOB_ID  $SGE_TASK_ID"
echo "=========================================================="


# Use an environment variable
n=$nval

echo "Running with n = $n"

python equilibrium_correlations.py

echo "Completed!"
