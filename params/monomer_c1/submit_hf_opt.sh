#!/bin/bash
# C1 monomer (beta-D-glucoside, exo core): ORCA HF/6-31G* optimisation from the RDKit start geometry.
# Runs in params/monomer_c1/hf_opt/ (git-ignored, see its .gitignore) so git pull/add cannot overwrite a running or finished job.
# On ARC:  cd /work/wrz135/callmann-tyler-md/params/monomer_c1 && bash submit_hf_opt.sh
# One core, 12 h limit. Expect roughly 2-4 h (the 3-O monomer took 14 cycles, 2 h 11 min, from a DFT start).
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p hf_opt
if grep -qs "ORCA TERMINATED NORMALLY" hf_opt/monomer_hf_opt.out; then
  echo "STOP: hf_opt/monomer_hf_opt.out already shows a finished run. Move hf_opt aside before rerunning."; exit 1
fi
if [ -s hf_opt/monomer_hf_opt.jobid ] && squeue -h -j "$(cat hf_opt/monomer_hf_opt.jobid)" 2>/dev/null | grep -q .; then
  echo "STOP: job $(cat hf_opt/monomer_hf_opt.jobid) is still queued/running."; exit 1
fi
n=$(head -1 dft-optimized-monomer.xyz | awk '{print $1}')
{ echo '! HF 6-31G* Opt TightSCF'; echo; echo '* xyz 0 1'; tail -n +3 dft-optimized-monomer.xyz; echo '*'; } > hf_opt/monomer_hf_opt.inp
c=$(awk '/^\* xyz/{f=1;next} /^\*$/{f=0} f&&NF==4{k++} END{print k+0}' hf_opt/monomer_hf_opt.inp)
[ "$n" = 61 ] && [ "$c" = "$n" ] || { echo "ERROR: expected 61 atoms, header $n, ORCA input has $c"; exit 1; }
cd hf_opt
jid=$(sbatch --parsable ../run_monomer_hf_opt.sbatch)
echo "$jid" > monomer_hf_opt.jobid
echo "submitted c1 HF opt, job $jid ($c atoms), working dir $(pwd)"
squeue -j "$jid" || true
