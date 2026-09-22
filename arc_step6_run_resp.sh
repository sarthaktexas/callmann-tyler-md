#!/usr/bin/env bash
set -euo pipefail

echo "=== ARC RESP pipeline step 6: two-stage RESP fit ==="
date
pwd
hostname

if [ ! -f monomer_resp.esp ]; then
  echo "ERROR: monomer_resp.esp is missing. Run step 5 first."
  exit 1
fi

if [ ! -f monomer_hf_opt.xyz ]; then
  echo "ERROR: monomer_hf_opt.xyz is missing."
  exit 1
fi

module load miniconda/24.4.0
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate ambertools
echo "resp: $(command -v resp)"

python3 make_resp_inputs.py --xyz monomer_hf_opt.xyz

resp -O \
  -i resp_stage1.in \
  -o resp_stage1.out \
  -p resp_stage1.pch \
  -q resp_stage1.qin \
  -t resp_stage1.qout \
  -e monomer_resp.esp

python3 make_resp_inputs.py --xyz monomer_hf_opt.xyz --stage1-qout resp_stage1.qout

resp -O \
  -i resp_stage2.in \
  -o resp_stage2.out \
  -p resp_stage2.pch \
  -q resp_stage2.qin \
  -t resp_stage2.qout \
  -e monomer_resp.esp

python3 make_resp_inputs.py --xyz monomer_hf_opt.xyz --stage2-qout resp_stage2.qout

echo
echo "=== Charge sum ==="
tail -1 resp_charges_stage2.dat
awk '/^[[:space:]]*[0-9]+/{sum+=$3} END{printf "awk charge_sum %.10f\n", sum; if (sum > 0.0001 || sum < -0.0001) exit 2}' resp_charges_stage2.dat

echo
echo "=== Outputs ==="
ls -lh resp_stage1.in resp_stage1.out resp_stage1.qout resp_stage2.in resp_stage2.out resp_stage2.qout resp_charges_stage2.dat

echo
echo "DONE step 6."
