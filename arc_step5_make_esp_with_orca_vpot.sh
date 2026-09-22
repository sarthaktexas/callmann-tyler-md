#!/usr/bin/env bash
set -euo pipefail

echo "=== ARC RESP pipeline step 5: ORCA vpot ESP generation ==="
date
pwd
hostname

if hostname | grep -qi '^login'; then
  echo "ERROR: Run this from an allocated compute shell, not a login node."
  echo "Run:"
  echo "  srun -p compute3 -n 1 -t 02:00:00 --cpus-per-task=1 --pty bash"
  exit 1
fi

module purge
module load orca

if [ ! -f monomer_hf_density.densities ]; then
  echo "monomer_hf_density.densities is missing."
  echo "This means the density job must be rerun with KeepDens."
  cat > monomer_hf_density.inp <<'EOF'
! HF 6-31G* TightSCF KeepDens

%output
  Print[P_Density] 1
end

* xyzfile 0 1 monomer_hf_opt.xyz
EOF
  jid="$(sbatch --parsable run_monomer_hf_density.sbatch)"
  echo "$jid" > monomer_hf_density.jobid
  echo "Submitted replacement density job id: $jid"
  echo "Wait for this job to finish, then rerun this step."
  exit 2
fi

if [ ! -f monomer_hf_density.gbw ]; then
  echo "ERROR: monomer_hf_density.gbw missing."
  exit 1
fi

if [ ! -f monomer_hf_opt.xyz ]; then
  echo "ERROR: monomer_hf_opt.xyz missing."
  exit 1
fi

python3 mk_resp_grid_from_orca.py --xyz monomer_hf_opt.xyz --points monomer_resp.vpot.xyz

echo "Available ORCA density-container entries:"
/apps/orca/6.0.1/orca_vpot monomer_hf_density.densities || true

/apps/orca/6.0.1/orca_vpot \
  monomer_hf_density.gbw \
  monomer_hf_density.scfp \
  monomer_resp.vpot.xyz \
  monomer_resp.vpot.out \
  monomer_hf_density

python3 mk_resp_grid_from_orca.py \
  --xyz monomer_hf_opt.xyz \
  --vpot-out monomer_resp.vpot.out \
  --resp-esp monomer_resp.esp \
  --convert

ls -lh monomer_resp.vpot.xyz monomer_resp.vpot.out monomer_resp.esp
head -3 monomer_resp.esp

echo
echo "DONE step 5."
