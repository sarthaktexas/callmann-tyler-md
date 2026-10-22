#!/bin/bash
# Make small, git-friendly copies of finished box trajectories WITH water, and run the analysis on ARC.
# On ARC, in a compute-node shell (srun -p compute1 -n 1 -c 2 -t 02:00:00 --pty bash):
#     cd /work/wrz135/callmann-tyler-md && git pull && bash scripts/reduce_box.sh sys_1to1 sys_2to1 chainonly_1 ...
# For each box NAME (needs boxes_c1/md_NAME.tpr and md_NAME.xtc) it writes analysis_data/NAME/:
#   md_1ns.xtc     whole system incl. water, molecules made whole, every FULL_STRIDE ps (default 1000 = 100 frames)
#   md.tpr         the run input (topology, masses, charges) for use with md_1ns.xtc
#   solute_200ps.xtc + solute.gro   pyrene + chains only, every SOLUTE_STRIDE ps (default 200), for the dynamics
#   result.csv, result_summary.txt  output of scripts/analyze_box.py (coverage, Rg, minimum distance)
# The full 10 ps trajectory stays on ARC in boxes_c1/md_NAME.xtc. Do NOT `git add .`: boxes_c1 holds multi-GB xtc files.
# After it finishes:  git add .gitattributes analysis_data && git commit -m "reduced C1 trajectories" && git push
set -eo pipefail
REPO="${REPO:-/work/wrz135/callmann-tyler-md}"; FULL_STRIDE="${FULL_STRIDE:-1000}"; SOLUTE_STRIDE="${SOLUTE_STRIDE:-200}"
if [ -z "$GMX" ]; then
  module purge; module load miniconda/24.4.0 >/dev/null 2>&1 || true
  eval "$(conda shell.bash hook)" 2>/dev/null; conda activate /work/wrz135/envs/gmx-cuda; GMX=$(command -v gmx)
fi
PY=python3; $PY -c "import numpy" 2>/dev/null || PY=""   # analysis needs numpy; if missing, run analyze_box.py locally on the reduced files
cd "$REPO/boxes_c1"
for S in "$@"; do
  case $S in
    sys_*to1) N=${S#sys_}; N=${N%to1}; NSOL=$((26 + 556 * N)); PYR=1 ;;
    chainonly_*) N=1; NSOL=556; PYR=0 ;;
    sys_pyronly) N=0; NSOL=26; PYR=1 ;;
    *) echo "unknown system $S"; continue ;;
  esac
  T=md_$S.tpr; X=md_$S.xtc; [ -f "$T" ] && [ -f "$X" ] || { echo "SKIP $S: $T / $X missing"; continue; }
  O="$REPO/analysis_data/$S"; mkdir -p "$O"; W=$(mktemp -d "${TMPDIR:-/tmp}/reduce.XXXXXX")
  echo "== $S  (solute atoms $NSOL, chains $N)  xtc $(du -h $X | cut -f1)"
  { echo "[ Solute ]"; seq -s ' ' 1 $NSOL; } > $W/idx.ndx
  echo System | $GMX trjconv -s $T -f $X -pbc whole -dt $FULL_STRIDE -o $O/md_1ns.xtc >/dev/null 2>&1
  cp $T $O/md.tpr
  echo System | $GMX trjconv -s $T -f $X -pbc whole -dt $SOLUTE_STRIDE -o $W/whole.xtc >/dev/null 2>&1
  echo Solute | $GMX trjconv -s $T -f $W/whole.xtc -n $W/idx.ndx -o $O/solute_200ps.xtc >/dev/null 2>&1
  echo Solute | $GMX trjconv -s $T -f $W/whole.xtc -n $W/idx.ndx -dump 0 -o $O/solute.gro >/dev/null 2>&1
  ls -la $O | awk 'NR>1{printf "   %s %s\n",$5,$9}'
  if [ -n "$PY" ] && [ "$N" -gt 0 ]; then
    FLAGS=""; [ $PYR = 0 ] && FLAGS="--no-pyrene"
    GMX=$GMX $PY $REPO/scripts/analyze_box.py --solute-only --gro $O/solute.gro --xtc $O/solute_200ps.xtc --nchains $N $FLAGS --out $O/result 2>&1 | tail -6 || echo "   (analysis failed for $S)"
    rm -rf $O/result_work
  else echo "   (analysis skipped: no python with numpy found, or pyrene-only box)"; fi
  rm -rf $W
done
du -sh "$REPO/analysis_data"
