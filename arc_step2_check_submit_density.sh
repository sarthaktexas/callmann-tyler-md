#!/usr/bin/env bash
set -euo pipefail

echo "=== ARC RESP pipeline step 2: check opt + submit density job ==="
date
pwd
hostname

if hostname | grep -qi '^login'; then
  echo "ERROR: ARC modules should be checked from an allocated compute shell, not a login node."
  echo "Run:"
  echo "  srun -p compute1 -n 1 -t 02:00:00 --cpus-per-task=8 --pty bash"
  echo "Then, inside that shell:"
  echo "  cd $(pwd)"
  echo "  bash arc_step2_check_submit_density.sh"
  exit 1
fi

if [ -f monomer_hf_opt.jobid ]; then
  jid="$(cat monomer_hf_opt.jobid)"
  echo "Recorded optimization job id: $jid"
  squeue -j "$jid" || true
  sacct -j "$jid" --format=JobID,JobName%24,State,ExitCode,Elapsed,MaxRSS || true
fi

echo
echo "=== ORCA optimization completion checks ==="
if [ ! -f monomer_hf_opt.out ]; then
  echo "ERROR: monomer_hf_opt.out is not present yet."
  exit 1
fi

grep -E "THE OPTIMIZATION HAS CONVERGED|ORCA TERMINATED NORMALLY|SCF CONVERGED" monomer_hf_opt.out || true

if ! grep -q "THE OPTIMIZATION HAS CONVERGED" monomer_hf_opt.out; then
  echo "ERROR: ORCA optimization has not cleanly converged yet."
  echo "Tail of monomer_hf_opt.out:"
  tail -60 monomer_hf_opt.out
  exit 1
fi

if ! grep -q "ORCA TERMINATED NORMALLY" monomer_hf_opt.out; then
  echo "ERROR: ORCA did not terminate normally."
  tail -60 monomer_hf_opt.out
  exit 1
fi

if [ ! -f monomer_hf_opt.xyz ]; then
  echo "ERROR: monomer_hf_opt.xyz was not produced."
  exit 1
fi

echo "Optimization converged and monomer_hf_opt.xyz exists."

echo
echo "=== Testing ORCA module load ==="
module purge
module load orca
echo "ORCA path: $(command -v orca || true)"
orca --version || true

echo
echo "=== Submitting density job ==="
jid="$(sbatch --parsable run_monomer_hf_density.sbatch)"
echo "$jid" > monomer_hf_density.jobid
echo "Submitted monomer_hf_density job id: $jid"
squeue -j "$jid" || true

echo
echo "DONE step 2."
