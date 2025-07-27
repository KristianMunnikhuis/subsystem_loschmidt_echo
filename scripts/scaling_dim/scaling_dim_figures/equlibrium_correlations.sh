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


# Get the argument from the command line
n=$1

python equilibrium_correlations.py

echo "Completed!"
