# C1 monomer (beta-D glucoside at C1, exo core): what exists and how to run it

Everything below mirrors the first (3-O) study, rebuilt on the C1 parameters with the corrected pyrene RESP charges. Written 21 Oct 2026.

## What was built (all checked locally; no C1 production MD has been run yet)
| item | location | status |
|---|---|---|
| monomer parameters (HF/6-31G* opt, PySCF-ESP standard RESP, GAFF2, residues MON/MNF/MNL/CPH/CPT) | `params/monomer_c1/` | done, all gates passed (`postopt.log`) |
| DP=10 chains s1-s4 (same cis/trans patterns as the 3-O chains: CCCCTTTTC, CCCTTTTTC, TTCCCTTCC, CTTCCCTTT) | `chains/c1_dp10/` | built by `scripts/build_chain.py`, 10 checks each |
| boxes: pyrene + 1/2/3/4 chains (Packmol 12-50 A shell, TIP3P 12 A), 4 chain-only (15 A), pyrene-only (28 A) | `boxes_c1/` | built by `scripts/build_box.py`; `grompp` accepts all nine |
| umbrella sampling v2, patterns 1-4 | `umbrella/v2_c1/` | prep stages tested locally for patterns 1 and 4 (build, solvate, index, EM, grompp of every stage) |
| per-frame analysis (pyrene coverage, Rg, min distance) | `scripts/analyze_box.py` | reproduces the earlier hand-computed numbers |

## Run on ARC
```
cd /work/wrz135/callmann-tyler-md && git pull
# 1. nine box jobs (100 ns each, one GPU; the 4:1 box is the slowest)
cd boxes_c1 && bash submit_all.sh                      # or: bash submit_all.sh sys_1to1 sys_2to1
# 2. umbrella (about 115 GPU-h per pattern; start with pattern 1, decide on the rest after looking at it)
cd ../umbrella/v2_c1 && mkdir -p logs
P=1; PREP=$(sbatch --parsable --export=ALL,PATTERN=$P 01_prep.sh)
sbatch --export=ALL,PATTERN=$P,DIR=in  --dependency=afterok:$PREP 02_windows.sh
sbatch --export=ALL,PATTERN=$P,DIR=out --dependency=afterok:$PREP 02_windows.sh
```
`02_windows.sh` exits 0 even if a window hits its time limit: check that each `windows_{in,out}/` has 30 `us_??.gro` before analysing.

## Analyse a finished box (any machine with GROMACS)
```
python scripts/analyze_box.py --tpr boxes_c1/md_sys_2to1.tpr --xtc boxes_c1/md_sys_2to1.xtc --nchains 2 --out analysis/c1_2to1
python scripts/analyze_box.py --tpr boxes_c1/md_chainonly_1.tpr --xtc boxes_c1/md_chainonly_1.xtc --nchains 1 --no-pyrene --out analysis/c1_chainonly_1
```
Quote the starting state with every number: pyrene coverage = 1 - SASA(pyrene, in complex)/SASA(pyrene, alone), same frame, no contact cutoff.
Every box has one pyrene at about 5 mM, so binding is forced by concentration; one run per box is a hypothesis, not a result.

## Rebuild from scratch
```
python scripts/build_chain.py --params params/monomer_c1 --out chains/c1_dp10 --tag s1 --like-old <old chain pdb>   # or --pattern / --seed
python scripts/build_box.py --kind nto1 --n 2 --out boxes_c1                                                    # also --kind chainonly --n K, --kind pyronly
```
Needs AmberTools (tleap, sander, ambpdb), Packmol, ParmEd and numpy.

## Known limits
* The monomer charges come from a single gas-phase HF/6-31G* conformer (as for the 3-O monomer); intramolecular sugar H-bonds can bias them.
* The first minimisation of each chain stops at its cycle limit (maximum force 0.2-1 kcal/mol/A); GROMACS EM follows in every box.
* Umbrella v2 has not run on ARC; its README lists what was and was not tested.
