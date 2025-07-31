#!/bin/bash -l

#$ -P fheating
#$ -m ea
#$ -N Equilibrium_Correlations
#$ -j y
module load python3/3.12.4

echo "================"
echo "Start date: $(date)"
echo "Job name : $JOB_NAME"
echo "Job ID : $JOB_ID  $SGE_TASK_ID"
echo "=========================================================="

for n in 1 2 3 4 5 6 7 8 9 10 11 12; do
    python equilibrium_connected_correlations.py "$n"
done 


echo "Completed!"

