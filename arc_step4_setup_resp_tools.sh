#!/usr/bin/env bash
set -euo pipefail

echo "=== ARC RESP pipeline step 4: setup/check RESP tools ==="
date
pwd
hostname

if hostname | grep -qi '^login'; then
  echo "ERROR: Run this from an allocated compute shell, not a login node."
  echo "Run:"
  echo "  srun -p compute3 -n 1 -t 02:00:00 --cpus-per-task=1 --pty bash"
  echo "Then, inside that shell:"
  echo "  cd $(pwd)"
  echo "  bash arc_step4_setup_resp_tools.sh"
  exit 1
fi

echo
echo "=== AmberTools conda env ==="
module load miniconda/24.4.0
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate ambertools
echo "ambertools" > ambertools_env.txt

echo
echo "=== AmberTools commands ==="
command -v resp
command -v antechamber
command -v respgen
command -v espgen || true
resp -h 2>&1 | head -20 || true

echo
echo "=== Multiwfn status ==="
module avail multiwfn 2>&1 | tee module_avail_multiwfn_step4.txt || true
if module load multiwfn 2>/dev/null; then
  echo "Multiwfn module loaded."
  command -v Multiwfn || command -v Multiwfn_noGUI || command -v multiwfn || true
else
  echo "No loadable Multiwfn module found."
  echo "Next recommended action: install/download a prebuilt Linux noGUI Multiwfn binary into the project or home directory."
  echo "Do not build Multiwfn from source yet; report this status first."
fi

echo
echo "DONE step 4."
