#!/bin/bash
#SBATCH --job-name=c1_hf_opt
#SBATCH --partition=compute3   # same partition as the 3-O monomer ORCA runs (proven); compute1 also works if compute3 is busy
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=12:00:00
#SBATCH --output=monomer_hf_opt.%j.out
#SBATCH --error=monomer_hf_opt.%j.err

set -euo pipefail

module purge
module load orca

export OMP_NUM_THREADS=1
orcabin="$(command -v orca)"
"$orcabin" monomer_hf_opt.inp > monomer_hf_opt.out
