#!/bin/bash
# On ARC: cd /work/wrz135/callmann-tyler-md/boxes_many && bash submit_many.sh [name ...]   (default: all five)
# Five independent GPU jobs of up to 24 h each (restartable: resubmit the same name to resume from md_<name>.cpt).
cd "$(dirname "$0")"; mkdir -p logs
LIST="${@:-many_p8_c8 many_p8_c16 many_p8_c24 many_p8_c32 many_p8_c0}"
for s in $LIST; do echo -n "$s: "; sbatch --parsable -J $s run_many.sh $s; done
