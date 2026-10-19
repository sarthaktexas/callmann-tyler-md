#!/usr/bin/env bash
# C1 monomer: everything after the ORCA HF/6-31G* optimisation. Run from params/monomer_c1 after monomer_hf_opt.out/.xyz exist.
#   bash run_postopt.sh              (refuses an unconverged optimisation)
#   ALLOW_UNCONVERGED=1 bash run_postopt.sh   plumbing test only; output must NOT be used for MD
# Needs on PATH: python with rdkit + pyscf (local conda env), AmberTools (antechamber, parmchk2, respgen, resp, tleap).
# PySCF ESP, standard Amber two-stage RESP, dipole gate, GAFF2, residues, tleap load test. Steps and gates: docs/PROTOCOL_glycomonomer_parameterization.md
set -euo pipefail
cd "$(dirname "$0")"
SMILES='C=C[C@@H]1C[C@H](C=C)[C@H]2C(=O)N(CCCCO[C@@H]3O[C@H](CO)[C@@H](O)[C@H](O)[C@H]3O)C(=O)[C@H]21'
NAT=61

echo "== 0. optimisation finished?"
if [ "${ALLOW_UNCONVERGED:-0}" != 1 ]; then
  grep -q "THE OPTIMIZATION HAS CONVERGED" monomer_hf_opt.out || { echo "STOP: no 'THE OPTIMIZATION HAS CONVERGED' in monomer_hf_opt.out"; exit 1; }
  grep -q "ORCA TERMINATED NORMALLY"       monomer_hf_opt.out || { echo "STOP: ORCA did not terminate normally"; exit 1; }
else echo "WARNING: ALLOW_UNCONVERGED=1 (test run, geometry is not an HF minimum)"; fi
grep -E "FINAL SINGLE POINT ENERGY" monomer_hf_opt.out | tail -1 || true
echo "optimisation cycles: $(grep -c 'GEOMETRY OPTIMIZATION CYCLE' monomer_hf_opt.out)"

echo "== G2. optimised geometry still the intended molecule (beta-D, exo, same stereo)"
python prep_monomer_geometry.py check monomer_hf_opt.xyz --smiles "$SMILES"
python prep_monomer_geometry.py mdl   monomer_hf_opt.xyz --smiles "$SMILES" -o monomer.mol

echo "== 4a. ESP (PySCF HF/6-31G*, Cartesian d)"
python esp_pyscf.py esp monomer_hf_opt.xyz -o monomer_resp.esp

echo "== 4b. standard two-stage RESP (Amber respgen)"
python3 -c "print('\n'.join(['0.0']*$NAT))" > q0.dat
antechamber -i monomer.mol -fi mdl -o mon.ac -fo ac -c rc -cf q0.dat -nc 0 -at gaff2 -rn MON -pf y > ac.log 2>&1
respgen -i mon.ac -o r1.in -f resp1
respgen -i mon.ac -o r2.in -f resp2
resp -O -i r1.in -o s1.out -p s1.pch -t s1.qout -e monomer_resp.esp
resp -O -i r2.in -o s2.out -p s2.pch -t s2.qout -q s1.qout -e monomer_resp.esp
python3 - <<'PY'
import numpy as np
q = np.array([float(x) for l in open('s2.qout') for x in l.split()])
el = [l.split()[0] for l in open('monomer_hf_opt.xyz').read().splitlines()[2:63]]
with open('resp_charges_stage2.dat', 'w') as f:
    f.write('# idx element charge  (standard Amber two-stage RESP, PySCF ESP)\n')
    for i, (e, c) in enumerate(zip(el, q), 1): f.write(f'{i:5d} {e:2s} {c: .8f}\n')
    f.write(f'# charge_sum {q.sum(): .10f}\n')
np.savetxt('q_final.dat', q, fmt='%.8f')
print('G5 charge sum %.2e' % q.sum())
PY

echo "== G4. dipole gate"
python esp_pyscf.py dipole monomer_hf_opt.xyz s2.qout --raw

echo "== 5. GAFF2 typing with final charges"
antechamber -i monomer.mol -fi mdl -o monomer_gaff2.mol2 -fo mol2 -c rc -cf q_final.dat -nc 0 -at gaff2 -rn MON -pf y > ac2.log 2>&1
parmchk2 -i monomer_gaff2.mol2 -f mol2 -o monomer.frcmod -s gaff2
echo "parmchk2 ATTN lines: $(grep -c ATTN monomer.frcmod)"

echo "== 6. residues from the molecular graph"
python build_residues.py monomer.mol monomer_gaff2.mol2 .

echo "== tleap load test: CPH-MNF-MON-MNL-CPT"
H=$(awk '/@<TRIPOS>HEADTAIL/{getline; print $1; getline; print $1}' res_MON.mol2 | tr '\n' ' ')
echo "head/tail atoms: $H"
cat > leap_test.in <<'LEAP'
source leaprc.gaff2
loadamberparams monomer.frcmod
MON = loadmol2 res_MON.mol2
MNF = loadmol2 res_MNF.mol2
MNL = loadmol2 res_MNL.mol2
CPH = loadmol2 res_CPH.mol2
CPT = loadmol2 res_CPT.mol2
chain = sequence { CPH MNF MON MNL CPT }
check chain
charge chain
saveamberparm chain leap_test.prmtop leap_test.inpcrd
quit
LEAP
tleap -f leap_test.in > leap_test.log 2>&1 || true
grep -E "Total unperturbed charge|Total charge|Warnings|Errors|FATAL|Unit is OK|missing|Could not" leap_test.log | head
echo "DONE post-opt"
