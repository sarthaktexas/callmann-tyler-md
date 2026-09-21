#!/bin/bash
#SBATCH --job-name=monomer_hf_den
#SBATCH --partition=compute3
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --output=monomer_hf_density.%j.out
#SBATCH --error=monomer_hf_density.%j.err

set -euo pipefail

module purge
module load orca

export OMP_NUM_THREADS=1
orcabin="$(command -v orca)"
"$orcabin" monomer_hf_density.inp > monomer_hf_density.out
