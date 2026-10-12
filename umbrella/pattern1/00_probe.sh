#!/bin/bash
#SBATCH --job-name=probe
#SBATCH --partition=gpu1a100
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=00:30:00
#SBATCH --output=/work/wrz135/callmann-tyler-md/umbrella/pattern1/logs/%x.%j.out
#SBATCH --error=/work/wrz135/callmann-tyler-md/umbrella/pattern1/logs/%x.%j.err
# Diagnostic only: does the job start, can it see /work, which GROMACS module works, does mdrun see the GPU.
# NOTE: Slurm silently fails (no log at all) if the --output directory does not exist. Run `mkdir -p logs` first.
echo "start $(date)"; echo "pwd=$PWD"; echo "user=$USER"
ls -ld /work/wrz135 /work/wrz135/callmann-tyler-md/umbrella/pattern1 || echo "WORK PATH PROBLEM"
touch /work/wrz135/callmann-tyler-md/umbrella/pattern1/logs/.write_test && echo "write ok" && rm /work/wrz135/callmann-tyler-md/umbrella/pattern1/logs/.write_test
source /work/wrz135/callmann-tyler-md/umbrella/pattern1/arc_env.sh
$GMX --version | grep -E "GROMACS version|GPU support|SIMD"
mkdir -p probe && cd probe
$GMX grompp -f ../mdp/em.mdp -c ../sys_pull.gro -p ../sys_pull.top -n ../index.ndx -o em.tpr -maxwarn 1 2>&1 | tail -3
$GMX mdrun -deffnm em -ntmpi 1 -ntomp 8 -nb gpu 2>&1 | tail -6
echo "EM exit: ${PIPESTATUS[0]}"
sed 's/nsteps = 50000/nsteps = 5000/' ../mdp/nvt.mdp > nvt.mdp
$GMX grompp -f nvt.mdp -c em.gro -p ../sys_pull.top -n ../index.ndx -o nvt.tpr -maxwarn 1 2>&1 | tail -3
$GMX mdrun -deffnm nvt $MD 2>&1 | tail -12
echo "NVT exit: ${PIPESTATUS[0]}"
echo "end $(date)"
