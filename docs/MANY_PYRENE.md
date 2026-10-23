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
