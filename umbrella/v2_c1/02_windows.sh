#!/bin/bash
#SBATCH --job-name=usc1_win
#SBATCH --partition=gpu1a100
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=08:00:00
#SBATCH --array=0-29%8
#SBATCH --output=/work/wrz135/callmann-tyler-md/umbrella/v2_c1/logs/%x.%A_%a.out
#SBATCH --error=/work/wrz135/callmann-tyler-md/umbrella/v2_c1/logs/%x.%A_%a.err
# Usage (from umbrella/v2_c1):  sbatch --export=ALL,PATTERN=2,DIR=in --dependency=afterok:<prep_jobid> 02_windows.sh   (and DIR=out)
# One window per array task (30 windows; near windows 80 ns, about 3 h; far windows 12 ns, about 30 min). Restartable.
set -eo pipefail
: "${DIR:?set DIR=in or DIR=out}"
V2=/work/wrz135/callmann-tyler-md/umbrella/v2_c1
source $V2/arc_env.sh
W=$(printf "%02d" $SLURM_ARRAY_TASK_ID)
cd windows_$DIR
[ -f us_$W.tpr ] || { echo "ERROR: windows_$DIR/us_$W.tpr missing - run 01_prep.sh first"; exit 1; }
if [ -f us_$W.gro ]; then echo "window $DIR/$W already finished"; exit 0; fi
if [ -f us_$W.cpt ]; then CPI="-cpi us_$W.cpt"; else CPI=""; fi
$GMX mdrun -deffnm us_$W $MD $CPI -cpt 15 -maxh 7.5 -px us_${W}_pullx -pf us_${W}_pullf
[ -f us_$W.gro ] && echo "WINDOW $DIR/$W DONE (pattern $PATTERN)" || echo "WINDOW $DIR/$W hit -maxh; resubmit: sbatch --export=ALL,PATTERN=$PATTERN,DIR=$DIR --array=$SLURM_ARRAY_TASK_ID 02_windows.sh"
