#!/bin/bash -l

#$ -P fheating
#$ -m ea
#$ -N Adiabatic_Projectors
#$ -j y

  # List of L values

echo "================"
echo "Start date: $(date)"
echo "Job name : $JOB_NAME"
echo "Job ID : $JOB_ID  $SGE_TASK_ID"
echo "=========================================================="

python adiabatic_projectors.py

echo "Completed!"
