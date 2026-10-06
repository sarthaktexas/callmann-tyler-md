#!/bin/bash
#SBATCH --job-name=pyrene_hf_opt
#SBATCH --partition=compute3
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=04:00:00
#SBATCH --output=pyrene_hf_opt.%j.out
#SBATCH --error=pyrene_hf_opt.%j.err

set -euo pipefail

module purge
module load orca

export OMP_NUM_THREADS=1
orcabin="$(command -v orca)"
"$orcabin" pyrene_hf_opt.inp > pyrene_hf_opt.out
