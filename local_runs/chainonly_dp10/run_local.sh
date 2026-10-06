#!/bin/bash
# Serial local run of chain-only controls 1-4 (DP=10): EM -> NVT 100 ps -> NPT 500 ps -> MD 10 ns. Restartable.
cd "$(dirname "$0")"
GMX=/Users/sarthakmohanty/.claude-science/conda/envs/md-local/bin/gmx
export OMP_NUM_THREADS=6
for N in 1 2 3 4; do
  S=chainonly_$N
  [ -f md_$N.gro ] && { echo "skip $N (done)"; continue; }
  [ -f em_$N.gro ]  || { $GMX grompp -f em.mdp -c $S.gro -p $S.top -o em_$N.tpr -maxwarn 1 && $GMX mdrun -deffnm em_$N -ntmpi 1 -ntomp 6; }
  [ -f nvt_$N.gro ] || { $GMX grompp -f nvt.mdp -c em_$N.gro -p $S.top -o nvt_$N.tpr -maxwarn 1 && $GMX mdrun -deffnm nvt_$N -ntmpi 1 -ntomp 6; }
  [ -f npt_$N.gro ] || { $GMX grompp -f npt.mdp -c nvt_$N.gro -t nvt_$N.cpt -p $S.top -o npt_$N.tpr -maxwarn 1 && $GMX mdrun -deffnm npt_$N -ntmpi 1 -ntomp 6; }
  [ -f md_$N.tpr ]  || $GMX grompp -f md.mdp -c npt_$N.gro -t npt_$N.cpt -p $S.top -o md_$N.tpr -maxwarn 1
  if [ -f md_$N.cpt ]; then $GMX mdrun -deffnm md_$N -ntmpi 1 -ntomp 6 -cpi md_$N.cpt -cpt 15
  else $GMX mdrun -deffnm md_$N -ntmpi 1 -ntomp 6 -cpt 15; fi
  echo "DONE chainonly_$N $(date)"
done
echo ALL_DONE
