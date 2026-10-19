# Umbrella sampling v2 (patterns 1-4), revised after the pattern-1 pilot

> **19 Oct 2026 note.** `01_prep.sh` now builds from `boxes/sys_1to1_pyrcorr.top` (corrected pyrene charges). The chains (`chains/dp10/chain_sN_min.pdb`, built from the original 3-O monomer) still have an L-glucose sugar and a non-norbornene core stereochemistry; see `docs/PROTOCOL_glycomonomer_parameterization.md`. Results from this protocol describe that model, not the intended D-glucose compound.

## Why v2
The pilot (`../pattern1`) was not converged. Its well depth fell with simulation time (10 ns blocks: -30, -31, -21, -13, -13 kJ/mol; all-data -20.6;
30-50 ns -13.0 +- 1.8). Bound windows (0.65-1.2 nm) had only about 17-120 independent samples each in 48 ns, windows started from a fast pull-out
(0.0012 nm/ps) of a wrapped complex, and start vs end snapshots showed pyrene coverage dropping 0.63 -> 0.44 in windows 0.85-1.3 nm
(see `analysis/pmf/pattern1/`). v2 addresses sampling and start-structure bias; it does not change the reaction coordinate.

## Changes relative to the pilot
| | pilot | v2 |
|---|---|---|
| start | final frame of a 100 ns bound run (pattern 1 only) | unbound chain_sN_min + pyrene at 3.1 nm, identical recipe for N = 1-4 |
| chain relaxation | none | 20 ns with pyrene restrained at 3.1 nm (extended starting chain collapses) |
| window seeding | one fast pull-out (0.0012 nm/ps) | slow pull-in (0.0002 nm/ps) **and** 5 ns hold at 0.55 nm + slow pull-out: two seed sets, `in` and `out` |
| windows | 28: 0.65-3.0 nm | 30: 0.55-1.0 (0.05), 1.1-1.8 (0.1), 1.9-3.0 (0.1). New 0.55, 0.60 nm windows bracket the well bottom |
| length | 50 ns each, 2 ns discarded | near (<1.85 nm, 18 windows): 80 ns, **first 30 ns discarded**; far (12 windows): 12 ns, first 2 ns discarded |
| pull reference | `pull-group1-pbcatom` only | plus `pull-pbc-ref-prev-step-com = yes` (a fresh chain starts extended; grompp refuses otherwise) |
| pull output | 0.2 ps | 1 ps |

Convergence check: compare the `in` and `out` PMFs (`analysis/pmf/run_v2.py`) **against the null spread measured by `analysis/pmf/null_calibration.py`**.
With no hysteresis at all and a 0.3 ns correlation time in the near windows (8 synthetic pairs, true depth 13.8 kJ/mol): depth sd 2.4 kJ/mol, |depth_in - depth_out|
up to 7.0, mean |W_in - W_out| median 1.3 (max 3.9). So one in/out pair detects only large hysteresis (> ~4 kJ/mol). Pilot bound windows had correlation times of about
2.6 ns, so real scatter will be larger: expect roughly +/- 4 kJ/mol on a depth from one direction. If the directions disagree, or each direction's own halves disagree,
do not quote the PMF; extend the near windows, pool replicates, or use enhanced sampling (replica exchange between windows needs an MPI GROMACS build, which the
conda-forge CUDA env does not provide; a second collective variable, e.g. a chain-pyrene contact count, is the other option).

## Facts checked before writing these scripts
* The four chains share one topology (`system1` section of `boxes/chainonly_{1..4}.top` is identical). Cis/trans pattern lives only in the starting
  coordinates (`chains/dp10/chain_s{1..4}_min.pdb`, same atom order). One `sys_pull.top` recipe therefore serves all patterns, and the CORE group and
  `pbcatom 369` are the same.
* Cis/trans identity is held by the double-bond torsion, not by anything checked in this directory; confirm junction dihedrals stay put in a finished window
  before interpreting patterns as different.
* Locally tested with GROMACS 2025.2 (CPU): build, solvation, index, `grompp` for em/nvt/npt/relax/pull_in/hold/pull_out and window mdp (patterns 1 and 4 all stages,
  2 and 3 build + em), the pull-length planner, and `make_windows.py` with mock pull data. NOT tested: any production run, `trjconv -pbc whole` on a real
  relaxed system, the Slurm scripts on ARC, whether pulling reaches 0.55 nm without sterics (the prep script fails loudly if a window has no frame within 0.03 nm).

## Run order on ARC (per pattern; submit from umbrella/v2)
```
cd /work/wrz135/callmann-tyler-md && git pull
cd umbrella/v2 && mkdir -p logs
P=2                                                                   # 1, 2, 3 or 4
PREP=$(sbatch --parsable --export=ALL,PATTERN=$P 01_prep.sh)          # ~2 h (about 50 ns of MD) on gpu1a100
sbatch --export=ALL,PATTERN=$P,DIR=in  --dependency=afterok:$PREP 02_windows.sh
sbatch --export=ALL,PATTERN=$P,DIR=out --dependency=afterok:$PREP 02_windows.sh
```
`02_windows.sh` exits 0 even if a window hits `-maxh 7.5`; check that each `windows_{in,out}/` has 30 `us_??.gro` before analysing and resubmit missing ones with
`sbatch --export=ALL,PATTERN=$P,DIR=in --array=<n> 02_windows.sh`. A near window needs about 3 h at the pilot speed.
Then copy `pattern$P/windows_*/us_??_pullx.xvg`, `windows_*.txt` back and run `python analysis/pmf/run_v2.py umbrella/v2/pattern$P`.
Do not commit trajectories without deciding on Git LFS (see below); only the `pullx` files are needed for the PMF.

## Cost (estimate from the pilot's measured 657-693 ns/day per window)
* Per direction: 18 x 80 ns + 12 x 12 ns = 1,584 ns. Per pattern, both directions: 3,168 ns, about 4.8 GPU-days = 115 GPU-h.
* Prep (relax 20 + pull-in ~12 + hold 5 + pull-out ~13 ns): about 50 ns, about 2 GPU-h.
* Four patterns: about 470 GPU-h in total. Against the TACC plan in memory (about 200 SU for the old design), this is about 2.3x larger.
  Lonestar6 packing (3 windows per node) would put it near 470 SU, close to the 500 SU startup limit. These are estimates.

## Git LFS note
`.gitattributes` tracks `*.xtc` with LFS, but `git lfs` is not installed on this machine, so the repo holds 133-byte pointer stubs, not trajectories (found while
analysing the pilot). Copy trajectories from ARC directly (rsync) or install git-lfs before committing xtc files.
