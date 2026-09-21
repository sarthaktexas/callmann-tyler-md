#!/bin/bash
#SBATCH --job-name=monomer_hf_den
#SBATCH --partition=compute3
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=01:00:00
#SBATCH --mem=8G
#SBATCH --output=monomer_hf_density.%j.out
#SBATCH --error=monomer_hf_density.%j.err

set -euo pipefail

module purge
module load orca

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"
orca monomer_hf_density.inp > monomer_hf_density.out
