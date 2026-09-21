#!/usr/bin/env bash
set -euo pipefail

echo "=== ARC RESP pipeline step 1: environment + HF opt submit ==="
date
pwd
hostname

if hostname | grep -qi '^login'; then
  echo "ERROR: ARC modules should be checked from an allocated compute shell, not a login node."
  echo "Run:"
  echo "  srun -p compute3 -n 1 -t 02:00:00 --cpus-per-task=4 --pty bash"
  echo "Then, inside that shell:"
  echo "  cd $(pwd)"
  echo "  bash arc_step1_env_submit_opt.sh"
  exit 1
fi

echo
echo "=== Files in working directory ==="
ls -lah

echo
echo "=== ORCA modules ==="
module avail orca 2>&1 | tee module_avail_orca.txt

echo
echo "=== Multiwfn modules ==="
module avail multiwfn 2>&1 | tee module_avail_multiwfn.txt || true

echo
echo "=== Conda environments, if conda exists ==="
if command -v conda >/dev/null 2>&1; then
  conda env list | tee conda_envs.txt
else
  echo "conda not currently on PATH" | tee conda_envs.txt
fi

echo
echo "=== Finding geometry ==="
geom=""
for f in dft-optimized-monomer.xyz *.xyz; do
  if [ -f "$f" ]; then
    geom="$f"
    break
  fi
done

if [ -z "$geom" ]; then
  echo "ERROR: no .xyz geometry found in $(pwd)"
  exit 1
fi

echo "Using geometry: $geom"
head -5 "$geom"
natoms="$(head -1 "$geom" | awk '{print $1}')"
echo "Atom count from XYZ header: $natoms"

echo
echo "=== Validating ORCA optimization input ==="
if [ ! -f monomer_hf_opt.inp ] || grep -q 'SCF=Tight' monomer_hf_opt.inp; then
  {
    echo '! HF 6-31G* Opt TightSCF PAL4'
    echo
    echo '* xyz 0 1'
    tail -n +3 "$geom"
    echo '*'
  } > monomer_hf_opt.inp
fi
awk 'BEGIN{n=0;inxyz=0} /^\* xyz /{inxyz=1;next} /^\*$/{if(inxyz){inxyz=0}} inxyz && NF==4{n++} END{print "Coordinate lines in ORCA input:", n}' monomer_hf_opt.inp

echo
echo "=== Testing ORCA module load in this shell ==="
module purge
module load orca
echo "ORCA path: $(command -v orca || true)"
orca 2>&1 | sed -n '1,12p' || true

echo
echo "=== Submitting HF optimization ==="
jid="$(sbatch --parsable run_monomer_hf_opt.sbatch)"
echo "$jid" > monomer_hf_opt.jobid
echo "Submitted monomer_hf_opt job id: $jid"

echo
echo "=== Queue status ==="
squeue -j "$jid" || true

echo
echo "DONE step 1."
