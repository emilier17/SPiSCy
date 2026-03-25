#!/bin/bash
#SBATCH --job-name={rule}
#SBATCH --cpus-per-task={threads}
#SBATCH --mem={resources.mem_mb}M
#SBATCH --time={resources.runtime}:00
#SBATCH --output=results/logs/slurm/{rule}_{jobid}.out
#SBATCH --error=results/logs/slurm/{rule}_{jobid}.err

{exec_jobscript}
