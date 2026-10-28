# Round 2 runs (built Oct 28)

Why: the first many-pyrene set showed about 85% coverage whenever chains >= pyrene, so the chains-per-pyrene answer lies on the
pyrene-rich side (not yet sampled), and the polymer effect on pyrene clustering was confounded with dilution. All boxes are built with
`scripts/build_many.py` (random unbound start, corrected pyrene charges, GAFF2/TIP3P), pass every builder check and a local grompp
(em.mdp and md.mdp, 0 warnings). Protocol identical to the first set (100 ns; DispCorr still off, kept for comparability).

| group | box | composition | box edge (nm) | built [pyrene] / [chain] (mM) | atoms | matched to |
|---|---|---|---|---|---|---|
| loading | many_p16_c8 | 16 pyrene + 8 chains | 10.65 | 22.1 / 11.1 | 99,904 | c8 box |
| loading | many_p24_c8 | 24 pyrene + 8 chains | 10.65 | 33.2 / 11.1 | 99,410 | c8 box |
| loading | many_p32_c8 | 32 pyrene + 8 chains | 10.65 | 44.2 / 11.1 | 98,895 | c8 box |
| controls | many_p8_c0_box8 | 8 pyrene | 10.65 | 11.1 / – | 104,260 | c8 (same pyrene start, seed 3808) |
| controls | many_p8_c0_box32 | 8 pyrene | 14.10 | 4.8 / – | 249,712 | c32 (same pyrene start, seed 3832) |
| chains | many_p0_c8 | 8 chains | 10.65 | – / 11.1 | 100,763 | c8 box |
| chains | many_p0_c32 | 32 chains | 14.10 | – / 11.4 | 235,145 | c32 box |
| (pre-formed) | many_p0_c16 | 16 chains | 11.05 | – | 110,371 | c16 box; submitted by submit_preformed.sh |
| reps | many_p8_c8_r2, many_p8_c16_r2 | as c8 / c16 | | | | same coordinates, new velocities |
| reps | sys_1to1_r2, sys_2to1_r2 (boxes_c1) | as 1:1 / 2:1 | | | | same coordinates, new velocities |

Loading boxes use `--pyr-sep 16 --rpyr 40/42` (pyrene centres >= 16 A apart; min pyrene-pyrene atom distance 9.8-12.7 A at start) because
16+ pyrene cannot be 30 A apart in the packing sphere. Box edges were hit with `--box-edge` (tleap buffer tuned in two passes).
Build-time concentrations rise by about 15-20% after NPT (c8 went 11.0 -> 13.0 mM).

Submit on ARC:

    cd /work/wrz135/callmann-tyler-md && git pull
    bash boxes_many/submit_round2.sh                 # 13 jobs: loading (3), controls (2), chains (2), reps (4 incl. 2 in boxes_c1)
    bash boxes_many/submit_preformed.sh              # if not already queued: many_p0_c16, then many_p8_c16_pre (afterok)

GPU cost: ~100k-atom boxes about 9-10 h each; the two ~240k-atom boxes about 18-22 h (c32 ran 134 ns/day). If a job hits 24 h,
resubmit by name: `sbatch -J <name> run_many.sh <name>`.

Bring back: `rsync` each md_<name>.xtc/.log/.gro/.tpr into boxes_many/ (or boxes_c1/ for sys_*_r2); `python scripts/build_site.py`
then picks them up (they are registered in the dashboard). Mechanism analyses: `scripts/mechanism_analysis.py` + `scripts/mechanism_summary.py`.
Umbrella sampling (binding free energy) remains staged in umbrella/v2_c1 (README there; ~115 GPU-h per pattern) and is not in this batch.
