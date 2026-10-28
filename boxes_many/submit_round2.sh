#!/bin/bash
# Round 2 (Oct 28): loading series, concentration-matched controls, chain-only aggregation controls, replicates.
# On ARC:  cd /work/wrz135/callmann-tyler-md && git pull && bash boxes_many/submit_round2.sh [group ...]
#   groups: loading controls chains reps   (default: all four)
# Every job is run_many.sh / run_box.sh (EM -> 100 ps NVT with fresh random velocities -> 500 ps NPT -> 100 ns), 24 h, restartable:
# resubmit the same name to resume from its checkpoint. Replicates (_r2) start from the same coordinates (symlinked .gro/.top)
# and differ only by NVT velocities (gen-seed = -1). The 16-chain chain-only box is submitted by submit_preformed.sh, not here.
set -e
cd "$(dirname "$0")"; mkdir -p logs ../boxes_c1/logs
G=${@:-loading controls chains reps}
sub() { echo -n "$1: "; sbatch --parsable -J $1 run_many.sh $1; }
for g in $G; do case $g in
  loading)  for n in many_p16_c8 many_p24_c8 many_p32_c8; do sub $n; done ;;                       # 8 chains + 16/24/32 pyrene (c8 = 8 pyrene)
  controls) for n in many_p8_c0_box8 many_p8_c0_box32; do sub $n; done ;;                          # 8 pyrene alone in the c8 / c32 box volume
  chains)   for n in many_p0_c8 many_p0_c32; do sub $n; done ;;                                    # 8 / 32 chains, no pyrene (box of c8 / c32)
  reps)     for n in many_p8_c8_r2 many_p8_c16_r2; do sub $n; done
            ( cd ../boxes_c1 && for n in sys_1to1_r2 sys_2to1_r2; do echo -n "$n: "; sbatch --parsable -J c1_$n run_box.sh $n; done ) ;;
  *) echo "unknown group $g"; exit 1 ;;
esac; done
