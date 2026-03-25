#!/bin/bash
#SBATCH --job-name=spiscy_master_
#SBATCH --output=spiscy_master_.%j.out
#SBATCH --error=spiscy_master_.%j.err
#SBATCH --time=10:00:00 
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --account=your_account_name

# ============================
# SPiSCy job manager
# ============================

# This allows Snakemake to manage the pipeline jobs in a compute node (not the login node)

# Go to the workflow directory
cd /path to spiscy folder

# Run spiscy launcher script
./run_spiscy.sh