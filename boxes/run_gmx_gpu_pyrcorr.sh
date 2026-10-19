#!/bin/bash
#SBATCH --partition=gpu1a100
#SBATCH --chdir=/work/wrz135/callmann-tyler-md/boxes
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
#SBATCH --output=gmx_%x.%j.out
#SBATCH --error=gmx_%x.%j.err

# Same as run_gmx_gpu.sh, but uses sys_Nto1_pyrcorr.top (pyrene H RESP corrected, scripts/swap_pyrene_charges.py) and tags every output _pcN
# so nothing from the earlier runs is reused or overwritten. Start = the same Packmol box (sys_Nto1.gro); 100 ns production (md.mdp).
set -eo pipefail
N=${1:?usage: sbatch -J pc<N> run_gmx_gpu_pyrcorr.sh <1|2>}
cd /work/wrz135/callmann-tyler-md/boxes
SYS=sys_${N}to1

echo "box ${SYS} with corrected pyrene RESP charges (${SYS}_pyrcorr.top); outputs tagged _pc${N}"
echo "node: $(hostname)"
lscpu | grep "Model name"
nvidia-smi -L

GMX=""
# 1) conda-forge CUDA GROMACS (same env the umbrella pilot used); the site gromacs modules on ARC are not CUDA-usable
GMX_CONDA_ENV=/work/wrz135/envs/gmx-cuda
if [ -d "$GMX_CONDA_ENV" ]; then
  module purge
  module load miniconda/24.4.0 >/dev/null 2>&1 || true
  eval "$(conda shell.bash hook)" 2>/dev/null && conda activate "$GMX_CONDA_ENV" 2>/dev/null
  g=$(command -v gmx || true)
  if [ -n "$g" ] && "$g" --version 2>&1 | grep -qi "GPU support: *CUDA" \
     && "$g" grompp -f em.mdp -c ${SYS}.gro -p ${SYS}_pyrcorr.top -o probe_pc${N}.tpr -maxwarn 1 >/dev/null 2>&1; then
    GMX=$g; echo "using conda env $GMX_CONDA_ENV ($g)"
  else echo "conda env present but not usable (no gmx, no CUDA, or grompp failed)"; fi
fi
# 2) fallback: site modules
if [ -z "$GMX" ]; then
for m in gromacs/2025.2 gromacs/2023.3-gpu; do
  module purge
  module load $m >/dev/null 2>&1 || continue
  g=$(command -v gmx || command -v gmx_mpi || true)
  [ -n "$g" ] || continue
  if ldd "$g" 2>/dev/null | grep -q "not found"; then echo "skip $m: missing libs"; continue; fi
  if ! "$g" --version 2>&1 | grep -qi "GPU support: *CUDA"; then echo "skip $m: no CUDA"; continue; fi
  GMX=$g; echo "using $m ($g)"; break
done
fi
[ -n "$GMX" ] || { echo "ERROR: no usable CUDA GROMACS (conda env or module)"; exit 1; }

MD="-ntomp 8 -nb gpu -pme gpu -bonded gpu"

if [ ! -f em_pc${N}.gro ]; then
  $GMX grompp -f em.mdp  -c ${SYS}.gro -p ${SYS}_pyrcorr.top -o em_pc${N}.tpr -maxwarn 1
  $GMX mdrun -deffnm em_pc${N} -ntomp 8 -nb gpu
fi
if [ ! -f nvt_pc${N}.gro ]; then
  $GMX grompp -f nvt.mdp -c em_pc${N}.gro -p ${SYS}_pyrcorr.top -o nvt_pc${N}.tpr -maxwarn 1
  $GMX mdrun -deffnm nvt_pc${N} $MD
fi
if [ ! -f npt_pc${N}.gro ]; then
  $GMX grompp -f npt.mdp -c nvt_pc${N}.gro -t nvt_pc${N}.cpt -p ${SYS}_pyrcorr.top -o npt_pc${N}.tpr -maxwarn 1
  $GMX mdrun -deffnm npt_pc${N} $MD
fi
if [ ! -f md_pc${N}.tpr ]; then
  $GMX grompp -f md.mdp  -c npt_pc${N}.gro -t npt_pc${N}.cpt -p ${SYS}_pyrcorr.top -o md_pc${N}.tpr -maxwarn 1
fi
if [ -f md_pc${N}.cpt ]; then
  $GMX mdrun -deffnm md_pc${N} $MD -cpi md_pc${N}.cpt -cpt 15
else
  $GMX mdrun -deffnm md_pc${N} $MD -cpt 15
fi
echo "DONE ${N}to1 pyrcorr"
