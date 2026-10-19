#!/bin/bash
# C1 monomer (beta-D-glucoside, exo core): write the ORCA input from the RDKit start geometry and submit the HF/6-31G* optimisation.
# Run on ARC from a login or compute shell:  cd /work/wrz135/callmann-tyler-md/params/monomer_c1 && bash submit_hf_opt.sh
# The 3-O monomer needed 14 cycles / 2 h 11 min on 1 core from a DFT start; this start is MMFF, so allow more (limit 12 h).
set -euo pipefail
cd "$(dirname "$0")"
n=$(head -1 dft-optimized-monomer.xyz | awk '{print $1}')
{ echo '! HF 6-31G* Opt TightSCF'; echo; echo '* xyz 0 1'; tail -n +3 dft-optimized-monomer.xyz; echo '*'; } > monomer_hf_opt.inp
c=$(awk '/^\* xyz/{f=1;next} /^\*$/{f=0} f&&NF==4{k++} END{print k+0}' monomer_hf_opt.inp)
[ "$n" = 61 ] && [ "$c" = "$n" ] || { echo "ERROR: expected 61 atoms, header $n, ORCA input has $c"; exit 1; }
jid=$(sbatch --parsable run_monomer_hf_opt.sbatch)
echo "$jid" > monomer_hf_opt.jobid
echo "submitted c1 HF opt, job $jid ($c atoms)"
squeue -j "$jid" || true
