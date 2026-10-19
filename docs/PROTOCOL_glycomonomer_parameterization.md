# Protocol: parameterizing a ROMP glycomonomer (GAFF2 + HF/6-31G* RESP)

Worked example: **β-D-glucopyranoside attached at C1, exo norbornene-imide core, butyl spacer** ("C1 monomer").
Written 19 Oct 2026. Everything marked **[run]** was executed while writing this; everything marked **[not run]** was not.

Files used below: [prep_monomer_geometry.py](../params/monomer_c1/prep_monomer_geometry.py), [esp_pyscf.py](../params/monomer_c1/esp_pyscf.py), [c1_beta_D_exo_start.xyz](../params/monomer_c1/c1_beta_D_exo_start.xyz), [pyrene_resp_charges_corrected.dat](../params/pyrene/pyrene_resp_charges_corrected.dat).
Put the two scripts in `params/<compound>/` next to `mk_resp_grid_from_orca.py` (`esp_pyscf.py` imports its grid code).

---

## 0. Why this protocol exists: what went wrong in the first pass

All of these were found on 19 Oct 2026 from the files in the repo. Each has a gate in this protocol so it cannot recur silently.

| # | Finding | Evidence | Status |
|---|---|---|---|
| 1 | **Sugar is L-glucose** in both the old 3-O xyz and the new C1 xyz. | Chirality-aware match against β-D-glucopyranose (CIP 2R,3R,4S,5S,6R); both files match only the mirror image. | **[run]** |
| 2 | **Core stereochemistry is impossible for a norbornene-derived unit** in both files. C1 xyz: trans-fused imide (fusion H–C–C–H = 176°). Old xyz: cis-fused, but the two vinyls are trans to each other. | A real ring-opened norbornene imide is cis-fused with cis-1,3 vinyls (exo or endo). I embedded all norbornene stereoisomers; only two are unstrained, neither matches either file. | **[run]** |
| 3 | **The pyrene ESP file is ~1.96× too large.** | PySCF HF/6-31G* ESP at the same 1875 points: correlation 0.993, scale 1.956. PySCF check is sound (Tr(PS) = 106 = Σ Z, far-field V ≈ 1e-5). | **[run]** |
| 4 | Consequence of 3: pyrene H RESP charges were **+0.24 to +0.27** (ΣH = 2.55 e). With a correct ESP and the same fit they are **+0.14 to +0.17** (ΣH = 1.55 e); AM1-BCC is +0.133. | Refit with Amber `resp`. The custom H-symmetry function was not the cause: Amber's own `respgen` stage 2 on the old ESP gives the same inflated values. | **[run]** |
| 5 | The **monomer** ESP file is fine. | PySCF vs the file: correlation 0.999997, scale 0.9989. RESP-charge dipole 7.58 D vs QM 7.61 D. | **[run]** |
| 6 | The repo's custom stage 2 (`make_resp_inputs.py`) is **not the standard Amber two-stage RESP**. It re-fits every atom (ivary = 0); Amber's `respgen` freezes everything except CH2/CH3 groups at stage-1 values (ivary = −99). On the old monomer ESP the two protocols differ by rms 0.049 e, max 0.231 e (atom 27, a carbon bonded to two O: 0.236 vs 0.467). Both reproduce the QM dipole, so the dipole gate cannot tell them apart. | Same ESP, both protocols. | **[run]**. Which is "better" is not established; the standard one is the published method. |

Not established: why `orca_vpot` produced a bad ESP for pyrene but a good one for the monomer, using identical scripts. Do not assume the ORCA vpot step is safe — use G3 below every time.

What this implies for existing MD (not assessed, hypothesis only): all runs used the inflated pyrene charges; the old monomer used the L-sugar, trans-1,3 core and non-standard stage 2.

---

## 1. Define the molecule as a stereo-defined SMILES (ChemDraw and Rowan are optional)

The identity of the molecule is its SMILES. The 3D geometry is only a starting point (it is re-optimised at HF/6-31G*), so a Rowan DFT geometry is not required. ChemDraw is only one way to write the SMILES; you can type it.

Decisions to record for every new monomer:

1. Sugar: **D**-gluco (natural), pyranose.
2. Attachment: C1 (glycoside) or C3 (ether); anomer **α/β** for C1.
3. Core: **exo or endo** (ring-opened norbornene imide, cis-fused, cis-1,3-divinyl).
4. Spacer length and N-substituent.

**C1 monomer chosen here: β-D, exo.**

```
C=C[C@@H]1C[C@H](C=C)[C@H]2C(=O)N(CCCCO[C@@H]3O[C@H](CO)[C@@H](O)[C@H](O)[C@H]3O)C(=O)[C@H]21
```

Formula C21H31NO8, 61 atoms (the vinyl-capped ring-opened model used to derive the repeat unit). Valid core references (N-methyl stand-ins): exo `C=C[C@@H]1C[C@H](C=C)[C@H]2C(=O)N(C)C(=O)[C@H]21`, endo `C=C[C@@H]1C[C@H](C=C)[C@@H]2C(=O)N(C)C(=O)[C@@H]21`.

If you do draw in ChemDraw: **wedge/hash every stereocentre** (the two vinyl-bearing ring carbons, the fusion carbons, all sugar carbons). A perspective chair with plain bonds carries no stereo; the original drawing had plain bonds at the vinyl-bearing carbons and none on the sugar. Then use *Copy As → SMILES* and run gate G1 on it.

## 2. Build the starting geometry  **[run]**

Environment: RDKit ≥ 2022.09 (conda-forge).

```bash
python prep_monomer_geometry.py build "<SMILES>" -o start.xyz     # 300 ETKDG conformers, lowest MMFF94 energy
python prep_monomer_geometry.py check start.xyz --smiles "<SMILES>"
```

**Gate G1 (must pass, exit code 0):** formula C21H31NO8; sugar reported as `D-gluco, 1-O-alkyl, beta`; core reported as `exo`; stereo identical to the input SMILES. The same `check` command accepts any xyz (Rowan, ORCA, …) and will flag L-sugars and invalid cores — on the two original files it reports `L-gluco` and `INVALID` core.

Caveat: the lowest MMFF conformer is arbitrary for a flexible molecule, and RESP charges depend on conformation (glucose can fold into intramolecular H-bonds). For consistency with the old monomer we use a single conformer. If C1 and 3-O results differ, check conformer artefacts before concluding it is chemistry.

## 3. HF/6-31G* geometry optimisation on ARC (ORCA)

Use a new directory per compound (`params/monomer_c1/`, so nothing in `params/monomer/` is touched). The existing scripts hard-code `monomer_*` filenames, so copy them into the new directory unchanged and keep the names. Copy `arc_step1_env_submit_opt.sh`, `arc_step2_check_submit_density.sh`, the two `run_monomer_hf_*.sbatch` files and `mk_resp_grid_from_orca.py` unchanged (run conventions: `docs/ARC_RESP_RUN_SEQUENCE.md`).

```bash
# on ARC, compute-node shell (modules only work there)
srun -p compute3 -n 1 -t 02:00:00 --cpus-per-task=1 --pty bash
cd /work/wrz135/callmann-tyler-md/params/monomer_c1
cp start.xyz dft-optimized-monomer.xyz        # name expected by step 1
bash arc_step1_env_submit_opt.sh              # writes monomer_hf_opt.inp: ! HF 6-31G* Opt TightSCF, charge 0 mult 1
```

Reference: the old 61-atom monomer took 14 optimisation cycles and **2 h 11 min on 1 core**. The sbatch limit is 4 h. `squeue`/`sacct` then check `THE OPTIMIZATION HAS CONVERGED` and `ORCA TERMINATED NORMALLY` (step 2 does this check for you).

**Gate G2:** copy `monomer_hf_opt.xyz` back, then
```bash
python prep_monomer_geometry.py mdl monomer_hf_opt.xyz --smiles "<SMILES>" -o monomer.mol
```
This fails if the optimised geometry no longer matches the intended stereo/connectivity, and writes an MDL file with correct bond orders **in the same atom order** (needed so charges line up).

## 4. ESP and RESP charges

The ORCA density + `orca_vpot` steps (2, 5) are optional; step 1 (optimisation) is all you need from ORCA. Compute the ESP with PySCF instead:

```bash
pip install pyscf                                   # [run locally: pyrene ESP in 8 s; a 61-atom monomer SCF + ESP took ~40 s in my check]
python esp_pyscf.py esp monomer_hf_opt.xyz -o monomer_resp.esp
```
HF/6-31G* with Cartesian d, Merz-Kollman-style shells (scales 1.4/1.6/1.8/2.0 × vdW, 1.4 pt/Å²), Amber `.esp` format. The script asserts electron count = Σ Z and far-field V ≈ 0.
*Running PySCF on ARC [not run]:* probably `pip install pyscf` in the ambertools env on a compute node; the Mac route above is tested.

**Gate G3 (ESP scale).** If an ESP came from anywhere else (e.g. ORCA vpot), compare it with the PySCF ESP at the same points. Accept scale 0.98–1.02 and correlation > 0.999. The bad pyrene file gave 1.956 / 0.993.

Standard Amber two-stage RESP **[run]** (AmberTools; use this, not `make_resp_inputs.py`):

```bash
python3 -c "print('\n'.join(['0.0']*61))" > q0.dat          # placeholder charges, 61 = number of atoms
antechamber -i monomer.mol -fi mdl -o mon.ac -fo ac -c rc -cf q0.dat -nc 0 -at gaff2 -rn MON -pf y
respgen -i mon.ac -o r1.in -f resp1
respgen -i mon.ac -o r2.in -f resp2
resp -O -i r1.in -o s1.out -p s1.pch -t s1.qout -e monomer_resp.esp
resp -O -i r2.in -o s2.out -p s2.pch -t s2.qout -q s1.qout -e monomer_resp.esp   # final charges: s2.qout
```
`respgen` equivalences the H on each CH2/CH3 and freezes all other atoms at stage-1 values (ivary −99) — the Bayly et al. procedure.

**Gate G4 (dipole):** `python esp_pyscf.py dipole monomer_hf_opt.xyz s2.qout --raw` must report ratio 0.8–1.25 (old monomer: 1.00 with the repo charges, 1.01 with standard RESP). This catches a grossly scaled ESP **only for molecules with a dipole**; pyrene has none, which is why it slipped through. For non-polar molecules use G3.

**Gate G5:** Σq within 1e-4 of 0 (the corrected pyrene fit gives 2e-6); atoms that are symmetry-equivalent get equal charges; hydrogen charges should fall in the ranges seen for the old monomer under standard RESP **[run]**: sp3 C–H 0.03–0.14, sp2 C–H 0.13–0.19, O–H 0.37–0.47 (pyrene aromatic H after correction: 0.14–0.17). Values far outside, like the old pyrene H at 0.24–0.27, need investigation.

**Gate G6 (cross-check, [not run] here):** compare per-atom charges with AM1-BCC (`antechamber -c bcc`). Differences up to ~0.1 e on heavy atoms are normal; hydrogens differing by a factor > 1.5 should be investigated, as pyrene H did.

Corrected pyrene charges (same ESP method, standard RESP): [pyrene_resp_charges_corrected.dat](../params/pyrene/pyrene_resp_charges_corrected.dat), same atom order as `pyrene_hf_opt.xyz`. These have **not** been swapped into any topology.

## 5. GAFF2 typing  **[run on the C1 start geometry, zero placeholder charges]**

```bash
antechamber -i monomer.mol -fi mdl -o monomer_gaff2.mol2 -fo mol2 -c rc -cf <final charges one per line> -nc 0 -at gaff2 -rn MON -pf y
parmchk2 -i monomer_gaff2.mol2 -f mol2 -o monomer.frcmod -s gaff2
```
C1 monomer result: atom types c, c2, c3, c5, c6, h1, h2, ha, hc, ho, n, o, oh, os; parmchk2 flagged **0** `ATTN` (missing) parameters. `-cf` needs one charge per line in atom order.

## 6. Repeat-unit residues and chains  **[not done]**

`params/gaff2/make_residues.py` hard-codes the old monomer's atom names (`C5`, `H4`, `C9`, `H10`, `H11`, and the symmetric `pairs`). It would silently produce wrong residues on a different atom order. It must be rewritten to find the two vinyl cap carbons and the symmetric pairs from the molecular graph, then build MON / MNF / MNL / CPH / CPT. `scripts/junctions.py` loops over residues 2–10 and should work unchanged for DP = 10. Re-validate the chain with the same checks as before.

## 7. Order of work for the C1 monomer

1. Run G1 on the start geometry (done: passes).
2. ARC: step 1 optimisation (≈2–3 h serial) → G2.
3. Local: PySCF ESP → standard RESP → G3–G6 → GAFF2.
4. Write the graph-based residue builder; rebuild seed-1 chain (same cis/trans pattern as the 3-O chain).
5. Same protocol as the 3-O pilot: 1:1 box, chain-only control, umbrella (28 windows, k = 1000).
6. Before comparing C1 with 3-O: rebuild the 3-O monomer with the corrected stereo (D-sugar, valid core) and standard RESP, and use the corrected pyrene charges; otherwise the comparison mixes chemistry with errors.

## 8. Open items for the group

- Which stage-2 RESP convention to adopt for the project (recommended: Amber standard via `respgen`).
- Cause of the 1.96× ESP in the pyrene `orca_vpot` run (density file mismatch? stale file?). Not diagnosed.
- Whether to keep the vinyl-capped 61-atom model vs. a more chain-like capping for RESP; unchanged from the old pipeline.
- Anomer/core for any further monomers (α, endo, C3, other sugars): extend the reference table in `prep_monomer_geometry.py`.
