#!/bin/bash
#SBATCH --job-name=us_prep
#SBATCH --partition=gpu1a100
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=06:00:00
#SBATCH --output=/work/wrz135/callmann-tyler-md/umbrella/pattern1/logs/%x.%j.out
#SBATCH --error=/work/wrz135/callmann-tyler-md/umbrella/pattern1/logs/%x.%j.err
# EM -> 100 ps NVT -> 1 ns NPT -> 2 ns steered pull -> window start frames -> per-window tpr files.
set -eo pipefail
source /work/wrz135/callmann-tyler-md/umbrella/pattern1/arc_env.sh
G="-maxwarn 1"
[ -f em.gro ]  || { $GMX grompp -f mdp/em.mdp  -c sys_pull.gro -p sys_pull.top -n index.ndx -o em.tpr $G;  $GMX mdrun -deffnm em -ntmpi 1 -ntomp 8 -nb gpu; }
[ -f nvt.gro ] || { $GMX grompp -f mdp/nvt.mdp -c em.gro  -p sys_pull.top -n index.ndx -o nvt.tpr $G; $GMX mdrun -deffnm nvt $MD; }
[ -f npt.gro ] || { $GMX grompp -f mdp/npt.mdp -c nvt.gro -t nvt.cpt -p sys_pull.top -n index.ndx -o npt.tpr $G; $GMX mdrun -deffnm npt $MD; }
if [ ! -f pull.gro ]; then
  $GMX grompp -f mdp/pull.mdp -c npt.gro -t npt.cpt -p sys_pull.top -n index.ndx -o pull.tpr $G
  $GMX mdrun -deffnm pull $MD -px pull_pullx -pf pull_pullf
fi
python3 make_windows.py        # writes windows/conf_XX.gro, windows/umbrella_XX.mdp, windows.txt
for i in $(seq -f "%02g" 0 $(( $(wc -l < windows.txt) - 1 ))); do
  [ -f windows/us_$i.tpr ] || $GMX grompp -f windows/umbrella_$i.mdp -c windows/conf_$i.gro -p sys_pull.top -n index.ndx -o windows/us_$i.tpr -maxwarn 1
done
echo "PREP DONE: $(wc -l < windows.txt) windows"
