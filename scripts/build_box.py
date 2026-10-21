#!/usr/bin/env python3
"""Build a GROMACS-ready simulation box from the tleap residues, the same recipe as the first (3-O) box set.

  build_box.py --kind nto1 --n 2 --out boxes_c1                 pyrene fixed at the origin + chains s1..sN in a 12-50 A shell
  build_box.py --kind chainonly --n 1 --out boxes_c1            chain sN alone, 15 A TIP3P buffer
  build_box.py --kind pyronly --out boxes_c1                    pyrene alone, 28 A TIP3P buffer

Packmol (tolerance 2.0; chains inside sphere 50 A and outside sphere 12 A) -> TER between molecules -> tleap (GAFF2, TIP3P,
solvatebox 12 A) -> ParmEd to GROMACS .top/.gro. Pyrene uses params/pyrene/pyrene_gaff2_pyrcorr.mol2 (corrected RESP charges).
Needs packmol, tleap, parmed (AmberTools) on PATH. Prints a report; exit 1 if a check fails.
"""
import argparse, os, re, subprocess, sys

ap = argparse.ArgumentParser()
ap.add_argument('--kind', choices=['nto1', 'chainonly', 'pyronly'], required=True)
ap.add_argument('--n', type=int, default=1, help='chains per pyrene (nto1) or chain number (chainonly)')
ap.add_argument('--out', default='boxes_c1'); ap.add_argument('--chains', default='chains/c1_dp10'); ap.add_argument('--params', default='params/monomer_c1')
ap.add_argument('--pyrene-mol2', default='params/pyrene/pyrene_gaff2_pyrcorr.mol2'); ap.add_argument('--pyrene-pdb', default='params/pyrene/pyrene_gaff2.pdb')
ap.add_argument('--pyrene-frcmod', default='params/pyrene/pyrene.frcmod'); ap.add_argument('--packmol-seed', type=int, default=None)
a = ap.parse_args()
A = os.path.abspath
OUT, CH, P = A(a.out), A(a.chains), A(a.params); os.makedirs(OUT, exist_ok=True)
PM, PP, PF = A(a.pyrene_mol2), A(a.pyrene_pdb), A(a.pyrene_frcmod)
name = {'nto1': f'sys_{a.n}to1', 'chainonly': f'chainonly_{a.n}', 'pyronly': 'sys_pyronly'}[a.kind]
sh = lambda cmd: subprocess.run(cmd, shell=True, cwd=OUT, capture_output=True, text=True)
natoms = lambda p: sum(1 for l in open(p) if l.startswith(('ATOM', 'HETATM')))
ok = True
def check(c, m):
    global ok; print(('PASS ' if c else 'FAIL ') + m); ok &= bool(c)

res_loads = ''.join(f'{r} = loadmol2 {P}/res_{r}.mol2\n' for r in ('MON', 'MNF', 'MNL', 'CPH', 'CPT'))
chain_ff = f'loadamberparams {P}/monomer.frcmod\n{res_loads}'
pyr_ff = f'loadamberparams {PF}\nPYR = loadmol2 {PM}\n'
head = 'source leaprc.gaff2\nsource leaprc.water.tip3p\n'

if a.kind == 'nto1':
    chains = [f'{CH}/chain_s{k}_min.pdb' for k in range(1, a.n + 1)]
    sizes = [natoms(PP)] + [natoms(c) for c in chains]
    seed = a.packmol_seed or 1000 + a.n
    inp = f'tolerance 2.0\nfiletype pdb\nseed {seed}\noutput {name}_raw.pdb\n\nstructure {PP}\n  number 1\n  fixed 0. 0. 0. 0. 0. 0.\nend structure\n'
    for c in chains: inp += f'\nstructure {c}\n  number 1\n  inside sphere 0. 0. 0. 50.\n  outside sphere 0. 0. 0. 12.\nend structure\n'
    open(f'{OUT}/pack_{a.n}to1.inp', 'w').write(inp)
    r = sh(f'packmol < pack_{a.n}to1.inp > pack_{a.n}to1.log 2>&1')
    check(os.path.exists(f'{OUT}/{name}_raw.pdb') and 'Success!' in open(f'{OUT}/pack_{a.n}to1.log').read(), 'packmol finished successfully')
    atoms = [l for l in open(f'{OUT}/{name}_raw.pdb') if l.startswith(('ATOM', 'HETATM'))]
    check(len(atoms) == sum(sizes), f'raw pdb has {len(atoms)} atoms == {sum(sizes)}')
    out = []; k = 0
    for s in sizes: out += atoms[k:k + s] + ['TER\n']; k += s
    open(f'{OUT}/{name}.pdb', 'w').write(''.join(out) + 'END\n')
    leap = head + chain_ff + pyr_ff + f'sys = loadpdb {name}.pdb\ncheck sys\ncharge sys\nsolvatebox sys TIP3PBOX 12.0\nsaveamberparm sys {name}.prmtop {name}.inpcrd\nquit\n'
elif a.kind == 'chainonly':
    leap = head + chain_ff + f'sys = loadpdb {CH}/chain_s{a.n}_min.pdb\ncheck sys\ncharge sys\nsolvatebox sys TIP3PBOX 15.0\nsaveamberparm sys {name}.prmtop {name}.inpcrd\nquit\n'
else:
    leap = head + pyr_ff + f'sys = loadpdb {PP}\ncheck sys\ncharge sys\nsolvatebox sys TIP3PBOX 28.0\nsaveamberparm sys {name}.prmtop {name}.inpcrd\nquit\n'
open(f'{OUT}/solv_{name}.in', 'w').write(leap)
sh(f'tleap -f solv_{name}.in > leap_{name}.txt 2>&1')
log = open(f'{OUT}/leap_{name}.txt').read()
check('Errors = 0' in log, 'tleap errors = 0')
m = re.search(r'Total unperturbed charge:\s+(\S+)', log); check(m and abs(float(m.group(1))) < 1e-3, f'tleap net charge {m.group(1) if m else "?"}')
ncc = len(re.findall(r'Close contact', log)); print(f'  tleap close-contact warnings: {ncc}')
check(os.path.exists(f'{OUT}/{name}.prmtop'), 'prmtop written')

import parmed as pmd
s = pmd.load_file(f'{OUT}/{name}.prmtop', f'{OUT}/{name}.inpcrd')
s.save(f'{OUT}/{name}.top', overwrite=True); s.save(f'{OUT}/{name}.gro', overwrite=True)
top = open(f'{OUT}/{name}.top').read()
mol = re.search(r'\[ molecules \]\n(.*)', top, re.S).group(1)
print('  [ molecules ]:', ' | '.join(' '.join(l.split()) for l in mol.strip().split('\n') if not l.startswith(';')))
nwat = int(re.search(r'WAT\s+(\d+)', mol).group(1))
box = open(f'{OUT}/{name}.gro').read().split('\n')[-2].split()
print('  box (nm):', box, ' atoms:', len(s.atoms), ' waters:', nwat)
expected = {'nto1': 26 + 556 * a.n, 'chainonly': 556, 'pyronly': 26}[a.kind]
check(len(s.atoms) == expected + 3 * nwat, f'atom count {len(s.atoms)} == solute {expected} + 3 x {nwat} waters')
if a.kind != 'chainonly':
    q = [at.charge for at in s.atoms[:26]]
    check(abs(sum(q)) < 1e-4 and abs(max(q[16:]) - 0.17) < 0.01, f'pyrene charges corrected (max H {max(q[16:]):.3f}, sum {sum(q):.1e})')
print('RESULT:', 'ALL CHECKS PASSED' if ok else 'FAILED'); sys.exit(0 if ok else 1)
