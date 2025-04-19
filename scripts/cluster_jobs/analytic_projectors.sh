#!/bin/bash -l

#$ -P fheating
#$ -m ea
#$ -N Mutual_info_Analytic
#$ -j y

L_values=(18 20 22 24 26)  # List of L values

echo "================"
echo "Start date: $(date)"
echo "Job name : $JOB_NAME"
echo "Job ID : $JOB_ID  $SGE_TASK_ID"
echo "=========================================================="

for L in "${L_values[@]}"; do
    echo "Running with L = $L"
    julia analytic_projectors.jl $L
done

echo "Completed!"
