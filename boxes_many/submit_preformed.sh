#!/bin/bash
# On ARC: cd /work/wrz135/callmann-tyler-md/boxes_many && bash submit_preformed.sh
# Stage 1: 16-chain chain-only pre-assembly (100 ns). Stage 2 starts only if stage 1 exits OK (afterok).
# If stage 1 hits the 24 h limit, stage 2 stays pending (DependencyNeverSatisfied): scancel it, resubmit stage 1 by name
# (sbatch -J many_p0_c16 run_many.sh many_p0_c16, resumes from md_many_p0_c16.cpt), then sbatch -J many_p8_c16_pre run_preformed.sh.
cd "$(dirname "$0")"; mkdir -p logs
j1=$(sbatch --parsable -J many_p0_c16 run_many.sh many_p0_c16); echo "stage 1 many_p0_c16: $j1"
j2=$(sbatch --parsable --dependency=afterok:$j1 -J many_p8_c16_pre run_preformed.sh); echo "stage 2 many_p8_c16_pre: $j2 (after $j1)"
