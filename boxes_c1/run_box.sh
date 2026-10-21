#!/bin/bash
#SBATCH --partition=gpu1a100
#SBATCH --chdir=/work/wrz135/callmann-tyler-md/boxes_c1
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
#SBATCH --output=logs/%x.%j.out
#SBATCH --error=logs/%x.%j.err
# Usage (from boxes_c1, after mkdir -p logs):  sbatch -J c1_<name> run_box.sh <name>
#   name = sys_1to1 sys_2to1 sys_3to1 sys_4to1 chainonly_1..4 sys_pyronly      (or run submit_all.sh)
# EM -> 100 ps NVT -> 500 ps NPT -> 100 ns production (md.mdp). Restartable: finished stages are skipped, production resumes from md_<name>.cpt.
# Uses the conda CUDA GROMACS (the site gromacs modules are not CUDA-usable on ARC); -ntmpi 1 is required with -ntomp on that build.
set -eo pipefail
SYS=${1:?usage: sbatch -J c1_<name> run_box.sh <name>}
cd /work/wrz135/callmann-tyler-md/boxes_c1
[ -f $SYS.gro ] && [ -f $SYS.top ] || { echo "ERROR: $SYS.gro/.top not found in $(pwd)"; exit 1; }
echo "system $SYS   node: $(hostname)"; lscpu | grep "Model name"; nvidia-smi -L
GMX_CONDA_ENV=/work/wrz135/envs/gmx-cuda
module purge; module load miniconda/24.4.0 >/dev/null 2>&1 || true
eval "$(conda shell.bash hook)" 2>/dev/null && conda activate "$GMX_CONDA_ENV"
GMX=$(command -v gmx)
"$GMX" --version 2>&1 | grep -qi "GPU support: *CUDA" || { echo "ERROR: $GMX has no CUDA"; exit 1; }
"$GMX" --version 2>&1 | grep -E "GROMACS version|GPU support|SIMD instructions"
MD="-ntmpi 1 -ntomp 8 -nb gpu -pme gpu -bonded gpu"
[ -f em_$SYS.gro ]  || { $GMX grompp -f em.mdp  -c $SYS.gro -p $SYS.top -o em_$SYS.tpr -maxwarn 1;  $GMX mdrun -deffnm em_$SYS -ntmpi 1 -ntomp 8 -nb gpu; }
[ -f nvt_$SYS.gro ] || { $GMX grompp -f nvt.mdp -c em_$SYS.gro -p $SYS.top -o nvt_$SYS.tpr -maxwarn 1; $GMX mdrun -deffnm nvt_$SYS $MD; }
[ -f npt_$SYS.gro ] || { $GMX grompp -f npt.mdp -c nvt_$SYS.gro -t nvt_$SYS.cpt -p $SYS.top -o npt_$SYS.tpr -maxwarn 1; $GMX mdrun -deffnm npt_$SYS $MD; }
[ -f md_$SYS.tpr ]  || $GMX grompp -f md.mdp -c npt_$SYS.gro -t npt_$SYS.cpt -p $SYS.top -o md_$SYS.tpr -maxwarn 1
if [ -f md_$SYS.cpt ]; then $GMX mdrun -deffnm md_$SYS $MD -cpi md_$SYS.cpt -cpt 15; else $GMX mdrun -deffnm md_$SYS $MD -cpt 15; fi
echo "DONE $SYS"
