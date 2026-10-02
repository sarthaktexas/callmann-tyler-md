#!/bin/bash
#SBATCH --job-name=em_test
#SBATCH --partition=compute1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00
#SBATCH --output=em_test.%j.out
#SBATCH --error=em_test.%j.err

set -eo pipefail
module purge
module load gromacs/2025.2
cd /work/wrz135/callmann-tyler-md/boxes

echo "node: $(hostname)"
lscpu | grep -E "Model name" 
echo "avx512f: $(grep -c avx512f /proc/cpuinfo) cores report it"

gmx grompp -f em.mdp -c sys_1to1.gro -p sys_1to1.top -o em_1.tpr -maxwarn 1
gmx mdrun -deffnm em_1 -ntomp 4 -v
