#!/bin/bash
#SBATCH --job-name=us_win
#SBATCH --partition=gpu1a100
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=08:00:00
#SBATCH --array=0-27%8
#SBATCH --output=/work/wrz135/callmann-tyler-md/umbrella/pattern1/logs/%x.%A_%a.out
#SBATCH --error=/work/wrz135/callmann-tyler-md/umbrella/pattern1/logs/%x.%A_%a.err
# One window per array task (28 windows, at most 8 running at once). Restartable: resubmitting resumes from the checkpoint.
set -eo pipefail
source /work/wrz135/callmann-tyler-md/umbrella/pattern1/arc_env.sh
W=$(printf "%02d" $SLURM_ARRAY_TASK_ID)
cd windows
[ -f us_$W.tpr ] || { echo "ERROR: windows/us_$W.tpr missing - run 01_prep_pull.sh first"; exit 1; }
if [ -f us_$W.gro ]; then echo "window $W already finished"; exit 0; fi
if [ -f us_$W.cpt ]; then CPI="-cpi us_$W.cpt"; else CPI=""; fi
$GMX mdrun -deffnm us_$W $MD $CPI -cpt 15 -maxh 7.5 -px us_${W}_pullx -pf us_${W}_pullf
[ -f us_$W.gro ] && echo "WINDOW $W DONE" || echo "WINDOW $W hit -maxh; resubmit: sbatch --array=$SLURM_ARRAY_TASK_ID 02_umbrella.sh"
