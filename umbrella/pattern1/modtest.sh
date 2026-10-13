#!/bin/bash
#SBATCH --job-name=modtest
#SBATCH --partition=gpu1a100
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=00:25:00
#SBATCH --output=/work/wrz135/callmann-tyler-md/umbrella/pattern1/logs/modtest.%j.out
#SBATCH --error=/work/wrz135/callmann-tyler-md/umbrella/pattern1/logs/modtest.%j.err
P=/work/wrz135/callmann-tyler-md/umbrella/pattern1
mkdir -p $P/modtest && cd $P/modtest
hostname
lscpu | grep "Model name"
echo "avx512 flags on this CPU: $(grep -o 'avx512[a-z]*' /proc/cpuinfo | sort -u | tr '\n' ' ')"
nvidia-smi --query-gpu=name,driver_version --format=csv,noheader
for m in gromacs/2025.2 gromacs/2023.3-gpu gromacs/2023.3 gromacs/2022.4 gromacs/mpi/2022.4 gromacs/nonmpi/2022.4; do
  echo; echo "######## $m"
  module purge; module load $m 2>&1 | head -5
  g=$(command -v gmx || command -v gmx_mpi || true); echo "binary: ${g:-none}"
  [ -n "$g" ] || continue
  echo "ldd missing: $(ldd $g 2>/dev/null | grep 'not found' | tr -s ' ' | tr '\n' ';')"
  $g --version 2>&1 | grep -E "GROMACS version|MPI library|GPU support|SIMD instructions|CUDA (driver|runtime)|PLUMED"
  NT="-ntmpi 1"; case "$g" in *gmx_mpi) NT="";; esac
  $g grompp -f ../mdp/em.mdp -c ../sys_pull.gro -p ../sys_pull.top -n ../index.ndx -o t.tpr -maxwarn 1 > grompp.log 2>&1 \
    && echo "grompp ok" || { echo "grompp FAILED"; tail -4 grompp.log; continue; }
  for flag in "" "-nb gpu"; do
    timeout 240 $g mdrun -s t.tpr -deffnm t -nsteps 200 $NT -ntomp 8 $flag > mdrun.log 2>&1
    echo "mdrun [$flag] exit=$?"
    grep -E -i "illegal|fatal|cannot|no compatible|Performance" mdrun.log | head -3
  done
done
