#!/bin/bash
#SBATCH --job-name=usc1_prep
#SBATCH --partition=gpu1a100
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=08:00:00
#SBATCH --output=/work/wrz135/callmann-tyler-md/umbrella/v2_c1/logs/%x.%j.out
#SBATCH --error=/work/wrz135/callmann-tyler-md/umbrella/v2_c1/logs/%x.%j.err
# Usage (from umbrella/v2_c1, after mkdir -p logs):  sbatch --export=ALL,PATTERN=2 01_prep.sh
# build unbound start -> solvate -> EM -> NVT -> NPT -> 20 ns relax (pyrene held at 3.1 nm) -> slow pull-in -> 5 ns hold at 0.55 nm
# -> slow pull-out -> window start frames and tpr files for BOTH directions ("in" and "out").
# Topology source: boxes_c1/sys_1to1.top (C1 monomer, corrected pyrene charges); the chain topology is identical for all four cis/trans patterns.
# Restartable: every stage is skipped when its output exists. ~2 h on an A100 (about 50 ns of MD).
set -eo pipefail
V2=/work/wrz135/callmann-tyler-md/umbrella/v2_c1
source $V2/arc_env.sh        # cd's into $V2/pattern$PATTERN, sets $GMX and $MD
G="-maxwarn 1"
if [ ! -f sys_pull.gro ]; then
  python3 $V2/build_start.py $PATTERN . 7.0
  $GMX editconf -f solute_unboxed.gro -o solute_box.gro -box 7.0 7.0 7.0 -noc >/dev/null 2>&1
  $GMX solvate -cp solute_box.gro -cs spc216.gro -o solv_raw.gro -maxsol 100000 2>&1 | grep -E "Number of SOL|Generated" || true
  python3 - <<'PY'
lines=open('solv_raw.gro').read().split('\n'); n=int(lines[1]); body=lines[2:2+n]; box=lines[2+n]; new=[]
for l in body:
    if l[5:10].strip()=='SOL':
        nm={'OW':'O','HW1':'H1','HW2':'H2'}[l[10:15].strip()]; l=l[:5]+'WAT  '+f'{nm:>5}'+l[15:]
    new.append(l)
nw=sum(1 for l in new if l[5:10].strip()=='WAT' and l[10:15].strip()=='O')
open('sys_pull.gro','w').write('\n'.join([lines[0],lines[1]]+new+[box,''])); open('nwat.txt','w').write(str(nw)); print('waters',nw)
PY
  NW=$(cat nwat.txt)
  python3 - "$NW" <<'PY'
import sys; nw=sys.argv[1]; s=open('/work/wrz135/callmann-tyler-md/boxes_c1/sys_1to1.top').read(); i=s.index('[ molecules ]')
open('sys_pull.top','w').write(s[:i]+'[ molecules ]\n; Compound       #mols\nPYR                  1\nsystem1              1\nWAT               %s\n'%nw)
PY
  python3 $V2/make_index.py sys_pull.gro index.ndx
fi
M=$V2/mdp
[ -f em.gro ]  || { $GMX grompp -f $M/em.mdp  -c sys_pull.gro -p sys_pull.top -n index.ndx -o em.tpr $G;  $GMX mdrun -deffnm em -ntmpi 1 -ntomp 8 -nb gpu; }
[ -f nvt.gro ] || { $GMX grompp -f $M/nvt.mdp -c em.gro  -p sys_pull.top -n index.ndx -o nvt.tpr $G; $GMX mdrun -deffnm nvt $MD; }
[ -f npt.gro ] || { $GMX grompp -f $M/npt.mdp -c nvt.gro -t nvt.cpt -p sys_pull.top -n index.ndx -o npt.tpr $G; $GMX mdrun -deffnm npt $MD; }
[ -f relax.gro ] || { $GMX grompp -f $M/relax.mdp -c npt.gro -t npt.cpt -p sys_pull.top -n index.ndx -o relax.tpr $G; $GMX mdrun -deffnm relax $MD -px relax_pullx -pf relax_pullf; }
# stage helper: whole-molecule copy of a gro, then the pull length needed to reach a target distance
plan() { # <deffnm> <target nm> <rate>
  echo System | $GMX trjconv -s $1.tpr -f $1.gro -n index.ndx -pbc whole -o $1_whole.gro >/dev/null 2>&1
  python3 $V2/plan_pull.py $1_whole.gro index.ndx $2 $3
}
if [ ! -f pull_in.gro ]; then
  eval "$(plan relax 0.50 0.0002 | tr ' ' '\n' | sed 's/^/IN_/' | tr '\n' ' ')"; echo "pull-in: D0=$IN_D0 nsteps=$IN_NSTEPS"
  sed "s/@NSTEPS@/$IN_NSTEPS/" $M/pull_in.mdp > pull_in_run.mdp
  $GMX grompp -f pull_in_run.mdp -c relax.gro -t relax.cpt -p sys_pull.top -n index.ndx -o pull_in.tpr $G
  $GMX mdrun -deffnm pull_in $MD -px pull_in_pullx -pf pull_in_pullf
fi
[ -f hold.gro ] || { $GMX grompp -f $M/hold.mdp -c pull_in.gro -t pull_in.cpt -p sys_pull.top -n index.ndx -o hold.tpr $G; $GMX mdrun -deffnm hold $MD -px hold_pullx -pf hold_pullf; }
if [ ! -f pull_out.gro ]; then
  eval "$(plan hold 3.05 0.0002 | tr ' ' '\n' | sed 's/^/OUT_/' | tr '\n' ' ')"; echo "pull-out: D0=$OUT_D0 nsteps=$OUT_NSTEPS"
  sed "s/@NSTEPS@/$OUT_NSTEPS/" $M/pull_out.mdp > pull_out_run.mdp
  $GMX grompp -f pull_out_run.mdp -c hold.gro -t hold.cpt -p sys_pull.top -n index.ndx -o pull_out.tpr $G
  $GMX mdrun -deffnm pull_out $MD -px pull_out_pullx -pf pull_out_pullf
fi
for D in in out; do
  [ -f windows_$D.txt ] || python3 $V2/make_windows.py $D pull_$D
  for i in $(seq -f "%02g" 0 $(( $(wc -l < windows_$D.txt) - 1 ))); do
    [ -f windows_$D/us_$i.tpr ] || $GMX grompp -f windows_$D/umbrella_$i.mdp -c windows_$D/conf_$i.gro -p sys_pull.top -n index.ndx -o windows_$D/us_$i.tpr -maxwarn 1
  done
done
echo "PREP DONE: pattern $PATTERN, $(wc -l < windows_in.txt) windows in + $(wc -l < windows_out.txt) windows out"
