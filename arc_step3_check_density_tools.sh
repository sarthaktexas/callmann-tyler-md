#!/usr/bin/env bash
set -euo pipefail

echo "=== ARC RESP pipeline step 3: check density + tool availability ==="
date
pwd
hostname

if hostname | grep -qi '^login'; then
  echo "ERROR: ARC modules should be checked from an allocated compute shell, not a login node."
  echo "Run:"
  echo "  srun -p compute1 -n 1 -t 02:00:00 --cpus-per-task=8 --pty bash"
  echo "Then, inside that shell:"
  echo "  cd $(pwd)"
  echo "  bash arc_step3_check_density_tools.sh"
  exit 1
fi

if [ -f monomer_hf_density.jobid ]; then
  jid="$(cat monomer_hf_density.jobid)"
  echo "Recorded density job id: $jid"
  squeue -j "$jid" || true
  sacct -j "$jid" --format=JobID,JobName%24,State,ExitCode,Elapsed,MaxRSS || true
fi

echo
echo "=== ORCA density completion checks ==="
if [ ! -f monomer_hf_density.out ]; then
  echo "ERROR: monomer_hf_density.out is not present yet."
  exit 1
fi

grep -E "ORCA TERMINATED NORMALLY|SCF CONVERGED" monomer_hf_density.out || true

if ! grep -q "ORCA TERMINATED NORMALLY" monomer_hf_density.out; then
  echo "ERROR: ORCA density job did not terminate normally."
  tail -60 monomer_hf_density.out
  exit 1
fi

if [ ! -f monomer_hf_density.gbw ]; then
  echo "ERROR: monomer_hf_density.gbw was not produced."
  exit 1
fi

echo "Density job completed and monomer_hf_density.gbw exists."

echo
echo "=== ORCA helper tools ==="
module purge
module load orca
echo "orca_2mkl: $(command -v orca_2mkl || true)"
echo "orca_vpot: $(command -v orca_vpot || true)"
orca_2mkl monomer_hf_density -molden
ls -lh monomer_hf_density.gbw monomer_hf_density.molden.input

echo
echo "=== Multiwfn availability ==="
module avail multiwfn 2>&1 | tee module_avail_multiwfn_step3.txt || true
if module load multiwfn 2>/dev/null; then
  echo "Multiwfn module loaded."
  command -v Multiwfn || command -v multiwfn || true
else
  echo "No loadable Multiwfn module found."
fi

echo
echo "=== AmberTools availability ==="
if command -v conda >/dev/null 2>&1; then
  conda env list
fi
if command -v mamba >/dev/null 2>&1; then
  mamba env list || true
fi
echo "Current resp: $(command -v resp || true)"
echo "Current antechamber: $(command -v antechamber || true)"

echo
echo "DONE step 3."
