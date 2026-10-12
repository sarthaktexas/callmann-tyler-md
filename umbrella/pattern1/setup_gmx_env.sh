#!/bin/bash
# One-time install of a CUDA-enabled GROMACS from conda-forge, on the ARC login node (not a job).
# ARC's own gromacs modules were rejected by the probe job ("no CUDA" or missing libs; see README).
# Usage:  bash setup_gmx_env.sh list            # show the available CUDA builds first
#         bash setup_gmx_env.sh install '<build string>'   # e.g. the newest 2025.x nompi_cuda_12* build from the list
set -e
module load miniconda/24.4.0 2>/dev/null || true
ENV=/work/wrz135/envs/gmx-cuda
case "$1" in
  list)
    conda search -c conda-forge --override-channels "gromacs>=2025" 2>/dev/null | grep -i cuda | tail -15
    ;;
  install)
    B=${2:?give a build string from the list, e.g. gromacs=2025.2=nompi_cuda_12.6_hXXXX_0}
    mkdir -p /work/wrz135/envs
    conda create -y -p $ENV -c conda-forge --override-channels "$B"
    conda run -p $ENV gmx --version | grep -E "GROMACS version|GPU support|SIMD instructions|CUDA (driver|runtime)"
    ;;
  *) echo "usage: $0 list | install <build>"; exit 1;;
esac
