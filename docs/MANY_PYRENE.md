# Many-pyrene / many-chain boxes (boxes_many/)

Question: with several pyrene present, how many DP=10 C1 chains are needed per pyrene, and do pyrene stay dispersed, aggregate
with each other, or bind to chains? Fixed composition, random UNBOUND start, one replicate per box (hypothesis-grade).

| box | pyrene | chains | atoms | box (nm) | pyrene mM | chain mM | chain:pyrene |
|---|---|---|---|---|---|---|---|
| many_p8_c0  | 8 | 0  |  40,360 | 7.9  | 27.4 | 0    | pyrene-only control (aggregation reference) |
| many_p8_c8  | 8 | 8  | 100,686 | 10.6 | 11.0 | 11.0 | 1:1 |
| many_p8_c16 | 8 | 16 | 110,198 | 11.1 | 9.8  | 19.6 | 2:1 |
| many_p8_c24 | 8 | 24 | 211,024 | 13.6 | 5.3  | 16.0 | 3:1 |
| many_p8_c32 | 8 | 32 | 237,147 | 14.1 | 4.7  | 18.9 | 4:1 |

Chains k = 1..N use start coordinates chain_s{(k-1) mod 4 + 1}_min.pdb (cis/trans patterns s1-s4 equally; one shared topology).
Pyrene: params/pyrene/pyrene_gaff2_pyrcorr.mol2 (corrected RESP). Solute atoms come first: 8 x 26 pyrene, then N x 556 chain atoms.

## Build
    build_many.py --n-pyr 8 --n-chain N --out boxes_many        (c8, c16, c0: defaults R=50 A sphere)
    build_many.py --n-pyr 8 --n-chain 24 --radius 60 --rpyr 44  (R=50 did not pack: Packmol penalty stayed ~5e3)
    build_many.py --n-pyr 8 --n-chain 32 --radius 66 --rpyr 48
Packmol: pyrene fixed at random sites (centres >= 30 A apart, random orientation), chains >= 12 A (centre-to-atom) from every
pyrene centre. Measured start: closest pyrene-chain atom distance 7.5-7.6 A in every box; no pyrene-chain contact at t=0.
Packmol sphere (not a cube) because tleap solvatebox aligns the solute to its principal axes: a cubic 80 A region gave a 113 A
water box, a sphere gives a box that depends only on the radius.

## Caveats to state with any result
1. Concentration is NOT constant: the pyrene concentration falls 11.0 -> 4.7 mM from c8 to c32 because the packing sphere had to
   grow to fit more chains. Ratio and dilution are confounded. For a clean ratio scan, either hold the box fixed (needs a
   packing method that handles ~19 mM chains, e.g. pre-solvated random insertion) or add a different pyrene count per box.
2. Pyrene is on its own in c0 at 27 mM: expect aggregation; this is the reference for "pyrene-pyrene" contacts in the other boxes.
3. c24 and c32 are 211k / 237k atoms; a 100 ns run is likely close to the 24 h limit on one A100. Production is restartable;
   resubmit the same name. NSTEPS=<n> overrides md.mdp's 50,000,000 steps.
4. Single replicate, unbound start, 100 ns. State the averaging window with any number.
5. grompp (gmx 2025.2, em.mdp, -maxwarn 1 as for boxes_c1) passes locally for all five; md.mdp grompp passes for c8 and c32.
   Nothing has been run (EM/NVT/NPT/production) and nothing has been submitted to ARC.

## Run (ARC, same workflow as boxes_c1)
    cd /work/wrz135/callmann-tyler-md/boxes_many && bash submit_many.sh [name ...]

## Pre-formed micelle (order of addition): many_p0_c16 -> many_p8_c16_pre
Question: does pyrene added to an already-assembled polymer aggregate (as in the pyrene fluorescence assay) end up as buried
as pyrene that co-assembles with the chains (many_p8_c16)? Same composition (8 pyrene, 16 chains), same box, same force field.

Stage 1, many_p0_c16: 16 chains, no pyrene, random start, 100 ns (also the 16-chain crowding control the many set lacked).
    build_many.py --n-pyr 0 --n-chain 16 --radius 46 --buffer 8 --out boxes_many
  R and buffer were tuned to match c16's box: R=50/buffer 6 gave 11.70 nm and 41,560 waters (pyrene would be 8.3 mM);
  R=46 or 47 with buffer 6 gave 10.65 nm and 29.7k waters. R=46/buffer 8 gives 11.01 x 11.04 x 11.07 nm, 33,825 waters,
  against c16's 11.05 nm and 33,698 waters.
Stage 2, many_p8_c16_pre: scripts/load_pyrene.py takes the final frame md_many_p0_c16.gro and places 8 pyrene in the water,
  random orientation, every pyrene atom >= 1.5 nm from every chain atom, centres >= 3.0 nm apart (minimum image). Waters with
  any atom < 0.30 nm from a pyrene atom are removed (about 22 per pyrene). The output keeps c16's atom order (pyrene, chains,
  water), so it reuses many_p8_c16.top (only the WAT count changes) and analyze_many.py works unchanged. Then the usual
  EM -> 100 ps NVT (gen-vel = yes: fresh velocities) -> 500 ps NPT -> 100 ns via run_many.sh.
  The system1 moleculetype from the chain-only build matches c16's line for line except one improper written in reversed
  atom order (3-1-2-4 vs 4-2-1-3: the same dihedral). Shared atom types are identical.
Local tests (gmx 2025.2): grompp em.mdp and md.mdp pass for many_p0_c16; load_pyrene.py on the t=0 chain-only box (stand-in
  for the final frame) passes all checks (min pyrene-chain 1.54 nm, 173 waters removed, 9.9 mM), grompp em/md pass, and
  1000-step EM converges to Fmax < 1000 in 11 steps.
Run on ARC:
    cd /work/wrz135/callmann-tyler-md/boxes_many && git pull && bash submit_preformed.sh
  stage 2 waits on stage 1 (afterok). If stage 1 hits the 24 h limit: scancel stage 2, resubmit stage 1 by name, then
  sbatch -J many_p8_c16_pre run_preformed.sh.
Compare with many_p8_c16 over matched windows, stating that c16 starts from an unbound mixed state and _pre from an unbound
  state next to a pre-equilibrated aggregate. Extra readouts: pyrene depth from the aggregate COM, radial imide/glucose
  density (is it core-shell?), and whether pyrene enters the aggregate or stays at its surface. Single replicate each.
