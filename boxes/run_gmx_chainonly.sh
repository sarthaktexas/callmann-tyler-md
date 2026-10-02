#!/bin/bash
#SBATCH --partition=gpu1a100
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
#SBATCH --output=gmx_%x.%j.out
#SBATCH --error=gmx_%x.%j.err

set -eo pipefail
K=${1:?usage: sbatch run_gmx_chainonly.sh <1|2|3|4>}; N=co${K}
cd /work/wrz135/callmann-tyler-md/boxes
SYS=chainonly_${K}

echo "node: $(hostname)"
lscpu | grep "Model name"
nvidia-smi -L

GMX=""
for m in gromacs/2025.2 gromacs/2026.3 gromacs/2023.3-gpu gromacs/2024.2; do
  module purge
  module load $m >/dev/null 2>&1 || continue
  g=$(command -v gmx || command -v gmx_mpi || true)
  [ -n "$g" ] || continue
  if ldd "$g" 2>/dev/null | grep -q "not found"; then echo "skip $m: missing libs"; continue; fi
  if ! "$g" --version 2>&1 | grep -qi "GPU support: *CUDA"; then echo "skip $m: no CUDA"; continue; fi
  if ! "$g" grompp -f em.mdp -c ${SYS}.gro -p ${SYS}.top -o probe_${N}.tpr -maxwarn 1 >/dev/null 2>&1; then echo "skip $m: grompp failed"; continue; fi
  GMX=$g
  echo "using $m ($g)"
  break
done
[ -n "$GMX" ] || { echo "ERROR: no usable CUDA GROMACS module"; exit 1; }

MD="-ntomp 8 -nb gpu -pme gpu -bonded gpu"

if [ ! -f em_${N}.gro ]; then
  $GMX grompp -f em.mdp  -c ${SYS}.gro -p ${SYS}.top -o em_${N}.tpr -maxwarn 1
  $GMX mdrun -deffnm em_${N} -ntomp 8 -nb gpu
fi
if [ ! -f nvt_${N}.gro ]; then
  $GMX grompp -f nvt.mdp -c em_${N}.gro -p ${SYS}.top -o nvt_${N}.tpr -maxwarn 1
  $GMX mdrun -deffnm nvt_${N} $MD
fi
if [ ! -f npt_${N}.gro ]; then
  $GMX grompp -f npt.mdp -c nvt_${N}.gro -t nvt_${N}.cpt -p ${SYS}.top -o npt_${N}.tpr -maxwarn 1
  $GMX mdrun -deffnm npt_${N} $MD
fi
if [ ! -f md_${N}.tpr ]; then
  $GMX grompp -f md.mdp  -c npt_${N}.gro -t npt_${N}.cpt -p ${SYS}.top -o md_${N}.tpr -maxwarn 1
fi
if [ -f md_${N}.cpt ]; then
  $GMX mdrun -deffnm md_${N} $MD -cpi md_${N}.cpt -cpt 15
else
  $GMX mdrun -deffnm md_${N} $MD -cpt 15
fi
echo "DONE ${N}to1"
