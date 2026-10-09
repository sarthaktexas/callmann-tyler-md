# Sourced by the job scripts. Picks the first GROMACS module that is CUDA-enabled, has no missing libraries
# and can read the topology. Notes from earlier ARC runs: gromacs/2023.3 gmx_mpi links a missing PLUMED library;
# all builds are AVX-512 and crash on compute3 (AMD EPYC 7F32), so exclude that node set.
WORKDIR=/work/wrz135/callmann-tyler-md/umbrella/pattern1
cd "$WORKDIR"
echo "node: $(hostname)  job: ${SLURM_JOB_ID:-none}  array task: ${SLURM_ARRAY_TASK_ID:-none}"
lscpu | grep "Model name" || true
nvidia-smi -L || echo "WARNING: no nvidia-smi output"
GMX=""
for m in gromacs/2025.2 gromacs/2026.3 gromacs/2024.2 gromacs/2023.3-gpu; do
  module purge
  module load $m >/dev/null 2>&1 || { echo "skip $m: not loadable"; continue; }
  g=$(command -v gmx || command -v gmx_mpi || true)
  [ -n "$g" ] || { echo "skip $m: no gmx binary"; continue; }
  if ldd "$g" 2>/dev/null | grep -q "not found"; then echo "skip $m: missing libs"; continue; fi
  if ! "$g" --version 2>&1 | grep -qi "GPU support: *CUDA"; then echo "skip $m: no CUDA"; continue; fi
  GMX=$g; echo "using $m ($g)"; break
done
[ -n "$GMX" ] || { echo "ERROR: no usable CUDA GROMACS module"; exit 1; }
# -update gpu is left off because it is not used with COM pulling here; nb/pme/bonded on GPU
export GMX
MD="-ntmpi 1 -ntomp ${SLURM_CPUS_PER_TASK:-8} -nb gpu -pme gpu -bonded gpu"
