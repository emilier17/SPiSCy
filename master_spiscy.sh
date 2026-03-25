#!/bin/bash
#SBATCH --job-name=spiscy_master_SUR
#SBATCH --output=spiscy_master_SUR_.%j.out
#SBATCH --error=spiscy_master_SUR_.%j.err
#SBATCH --time=10:00:00 
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --account=def-pbegin

# ============================
# SPiSCy job manager
# ============================

# This allows Snakemake to manage the pipeline jobs in a compute node (not the login node)

# Go to the workflow directory
cd /home/emilier/projects/def-pbegin/emilier/spiscy_sur

# Run spiscy launcher script
./run_spiscy.sh