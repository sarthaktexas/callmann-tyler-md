#!/bin/bash
# On ARC: cd /work/wrz135/callmann-tyler-md/boxes_c1 && bash submit_all.sh [name ...]   (default: every C1 box)
# Nine independent GPU jobs of up to 24 h each; the 4:1 box (73k atoms) is the slowest.
cd "$(dirname "$0")"; mkdir -p logs
LIST="${@:-sys_1to1 sys_2to1 sys_3to1 sys_4to1 chainonly_1 chainonly_2 chainonly_3 chainonly_4 sys_pyronly}"
for s in $LIST; do echo -n "$s: "; sbatch --parsable -J c1_$s run_box.sh $s; done
