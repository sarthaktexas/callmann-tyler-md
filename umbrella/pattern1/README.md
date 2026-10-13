# Pyrene–chain umbrella sampling, pilot (chain pattern 1)

Goal: potential of mean force (PMF) for pulling one pyrene off one DP=10 glucose-norbornene-imide chain, giving
binding free energy ΔG_bind. Pilot checks that the reaction coordinate and window overlap work before the other three patterns.

## Starting system (already built, committed)
* `build_start.py`: final frame (100 ns) of the 1:1 seed-1 Modal run. Chain made whole, pyrene placed in the nearest image.
* `solvate.sh`: 7.0 nm cubic box, 11,022 TIP3P waters (WAT, O/H1/H2), 33,648 atoms. Topology `sys_pull.top` reuses the
  force-field section of `boxes/sys_1to1.top`; molecules: PYR 1, system1 1, WAT 11022. The system is neutral.
* `index.ndx`: extra groups `PYRENE` (atoms 1-26), `CHAIN` (27-582), `CORE` (120 core heavy atoms, `core_heavy_atoms.txt`), `SOLUTE`.
* Reaction coordinate: distance between COM of `CORE` and COM of `PYRENE` (pull geometry `distance`).
  `CORE` spans about 2.5 nm, so `pull-group1-pbcatom = 369` is set (GROMACS refuses otherwise).
  The chain is flexible, so this coordinate does not control chain conformation, which can be a source of hysteresis;
  the pilot is meant to show how large that problem is.
* Starting distance 0.655 nm. Windows: 0.65-1.00 nm every 0.05 nm, then 1.1-3.0 nm every 0.1 nm (28 windows), k = 1000 kJ/mol/nm², 50 ns per window (the ARC probe ran at about 660 ns/day on an A100, so a window takes roughly 2 h).
  Box size: at 3.0 nm the pyrene and the chain's far end stay more than 2 nm from their periodic images.

## What was tested locally (Mac, CPU, gromacs 2025.2)
* grompp for EM, NVT, pull and a window mdp: passes (one `-maxwarn 1`, same as the box scripts).
* Steepest-descent EM converges (Fmax < 1000 in 350 steps); 20 ps NVT and 20 ps pull at the real rate: no LINCS warnings.
* `make_windows.py` ran on the 20 ps test pull for two windows (0.65, 0.70 nm); the generated window mdp passes grompp.
* NOT tested: the NPT step, the 2 ns pull, windows above 0.73 nm, any window run, WHAM, and every slurm script on ARC.
  Local speed was too low (about 1 ns/day) to test further.
* ARC probe and prep jobs (Oct 9) started on gpu039/gpu040 (A100-PCIE-40GB, EPYC 7F32) and exited after 2-4 s: the module loop found
  no CUDA GROMACS. ARC has gromacs 2022.4, 2023.3, 2023.3-gpu, 2025.2, mpi/2022.4 and nonmpi/2022.4 (not 2024.2 or 2026.3). 2025.2 reported
  no CUDA; 2023.3-gpu reported missing libraries. Logs were written, so the earlier "missing logs directory" idea was wrong.
* Fixed: conda-forge CUDA GROMACS 2025.2 (`gromacs=2025.2=nompi_cuda_h7ac747b_1`, installed with CONDA_OVERRIDE_CUDA=12.0 from a compute node) in /work/wrz135/envs/gmx-cuda; `arc_env.sh` uses it first.
  Probe job 836069 on gpu040 (A100-PCIE-40GB, driver 13.20, EPYC 7F32): AVX2_256 build, EM exit 0, NVT exit 0, 657 ns/day for a 5,000-step run.

## One-time GROMACS install on ARC (login node)
```
cd /work/wrz135/callmann-tyler-md/umbrella/pattern1
bash setup_gmx_env.sh list                    # pick a 2025.x build with cuda_12 in its build string
bash setup_gmx_env.sh install '<build string>'
```

## Run order on ARC
```
cd /work/wrz135/callmann-tyler-md && git pull
cd umbrella/pattern1 && mkdir -p logs
sbatch 00_probe.sh                     # 30 min limit; read logs/probe.*.out / .err
sbatch 01_prep_pull.sh                 # EM, NVT, NPT, pull, window tpr files
sbatch --dependency=afterok:<prep_jobid> 02_umbrella.sh     # 28 windows as an array, 8 at a time
sbatch --dependency=afterok:<array_jobid> 03_wham.sh
```
Scripts assume the repo lives in `/work/wrz135/callmann-tyler-md` (taken from the earlier `boxes/` scripts). Edit `WORKDIR` in `arc_env.sh`
and the `--output` lines if that is wrong. A window that reaches `-maxh 7.5` is resumed by `sbatch --array=<n> 02_umbrella.sh`.

## Checks before believing the PMF
1. `wham/histo.xvg`: neighbouring windows should overlap. Add windows where they do not.
2. `wham/profile.xvg` should plateau beyond ~2.5 nm; ΔG_bind = plateau minus minimum.
3. `wham/bsres.xvg` gives bootstrap error bars; they assume uncorrelated samples, so treat them as lower bounds.
4. Compare the PMF from the first and second 10 ns of each window. A difference larger than a few kJ/mol means more sampling is needed.
5. The value is a free energy for one pyrene on one chain with a standard-state correction not applied. Do not read it as a chains-per-pyrene number.
