#!/bin/bash
#SBATCH --partition=gpu1a100
#SBATCH --chdir=/work/wrz135/callmann-tyler-md/boxes_many
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
#SBATCH --output=logs/%x.%j.out
#SBATCH --error=logs/%x.%j.err
# Pre-formed micelle experiment, stage 2: load 8 pyrene into the equilibrated 16-chain chain-only box, then run
# EM -> 100 ps NVT (fresh velocities, gen-vel = yes) -> 500 ps NPT -> 100 ns, exactly as run_many.sh.
#   Stage 1 (chain-only pre-assembly, also the 16-chain crowding control):  sbatch -J many_p0_c16 run_many.sh many_p0_c16
#   Stage 2 (this script):                                                   sbatch -J many_p8_c16_pre run_preformed.sh
#   or both at once with a dependency:                                       bash submit_preformed.sh
# Restartable: if many_p8_c16_pre.gro exists the loading step is skipped and run_many.sh resumes from its checkpoint.
# Comparator: many_p8_c16 (co-assembled, same composition, same PYR/system1/WAT moleculetypes, same atom order).
set -eo pipefail
CH=${CH:-many_p0_c16}; OUT=${OUT:-many_p8_c16_pre}; NCH=${NCH:-16}; NPYR=${NPYR:-8}; TMPL=${TMPL:-many_p8_c16}
cd /work/wrz135/callmann-tyler-md/boxes_many
if [ ! -f $OUT.gro ]; then
  [ -f md_$CH.gro ] || { echo "ERROR: md_$CH.gro not found: the chain-only run ($CH) has not finished"; exit 1; }
  grep -q "Finished mdrun" md_$CH.log || { echo "ERROR: md_$CH.log has no 'Finished mdrun' (resubmit $CH first)"; exit 1; }
  tail -n 60 md_$CH.log | grep -A1 -E "^ +Step +Time" | tail -1
  module purge; module load miniconda/24.4.0 >/dev/null 2>&1 || true
  python3 ../scripts/load_pyrene.py --gro md_$CH.gro --n-chain $NCH --n-pyr $NPYR \
      --template $TMPL.gro --top $TMPL.top --out $OUT | tee load_$OUT.log
  grep -q "ALL CHECKS PASSED" load_$OUT.log || { echo "ERROR: load_pyrene.py checks failed"; rm -f $OUT.gro $OUT.top; exit 1; }
fi
bash run_many.sh $OUT
