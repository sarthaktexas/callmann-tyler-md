#!/usr/bin/env python3
"""Build a many-pyrene / many-chain GROMACS box (random, unbound start) from the C1 tleap residues.

  build_many.py --n-pyr 8 --n-chain 16 --out boxes_many     8 pyrene + 16 DP=10 chains (2:1)
  build_many.py --n-pyr 8 --n-chain 0  --out boxes_many     pyrene-only control (same box)
  build_many.py --n-pyr 0 --n-chain 32 --out boxes_many     chain-only crowding control (same box)

Start state (unbound by construction):
  * pyrene: n_pyr molecules at random positions (min centre separation 30 A, all within RPYR of the origin), random orientation,
            fixed in Packmol
  * chains: random position/orientation inside a sphere of radius R (default 50 A, as in build_box.py), >= 12 A (centre-to-atom,
            same convention as build_box.py) from every pyrene centre; chain k uses chain_s{(k mod 4)+1}_min.pdb, so all four
            cis/trans patterns are used equally (the topology is identical for all patterns)
  * A sphere is used because tleap solvatebox aligns the solute to its principal axes: a cubic packing region of edge 80 A gave a
    113 A box, a sphere gives about 2R + 2*buffer for any arrangement, so the box is nearly the same for every composition.
  * Packmol tolerance 2.0 A; tleap GAFF2/TIP3P, solvatebox with an isotropic buffer; ParmEd -> GROMACS.
Solute atom order: all pyrene (26 each) first, then chains (556 each), then water.  Needs packmol, tleap, parmed on PATH.
"""
import argparse, math, os, random, re, subprocess, sys
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument('--n-pyr', type=int, required=True); ap.add_argument('--n-chain', type=int, required=True)
ap.add_argument('--radius', type=float, default=50.0, help='Packmol sphere radius (A); the TIP3P box is about 2*radius + 3 + 2*buffer')
ap.add_argument('--rpyr', type=float, default=36.0, help='pyrene centres lie within this radius (A)')
ap.add_argument('--buffer', type=float, default=6.0); ap.add_argument('--out', default='boxes_many')
ap.add_argument('--chains', default='chains/c1_dp10'); ap.add_argument('--params', default='params/monomer_c1')
ap.add_argument('--pyrene-mol2', default='params/pyrene/pyrene_gaff2_pyrcorr.mol2'); ap.add_argument('--pyrene-pdb', default='params/pyrene/pyrene_gaff2.pdb')
ap.add_argument('--pyrene-frcmod', default='params/pyrene/pyrene.frcmod'); ap.add_argument('--seed', type=int, default=None)
ap.add_argument('--name', default=None, help='output name (default many_p<NP>_c<NC>)')
ap.add_argument('--pyr-sep', type=float, default=30.0, help='minimum pyrene centre-centre distance (A); lower it for >8 pyrene')
ap.add_argument('--box-edge', type=float, default=None, help='target cubic box edge (nm): the solvatebox buffer is tuned (two tleap passes) to reach it')
a = ap.parse_args()
NP, NC, R = a.n_pyr, a.n_chain, a.radius
assert NP + NC > 0
A = os.path.abspath; OUT, CH, P = A(a.out), A(a.chains), A(a.params); os.makedirs(OUT, exist_ok=True)
PM, PP, PF = A(a.pyrene_mol2), A(a.pyrene_pdb), A(a.pyrene_frcmod)
name = a.name or f'many_p{NP}_c{NC}'
seed = a.seed or (3000 + 100 * NP + NC)
sh = lambda cmd: subprocess.run(cmd, shell=True, cwd=OUT, capture_output=True, text=True)
rec = lambda p: [l for l in open(p) if l.startswith(('ATOM', 'HETATM'))]
xyz = lambda l: np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])
ok = True
def check(c, m):
    global ok; print(('PASS ' if c else 'FAIL ') + m); ok &= bool(c)

# pyrene reference, recentred on its geometric centre so 'fixed' coordinates are unambiguous
pl = rec(PP); pc = np.mean([xyz(l) for l in pl], axis=0)
cen = []
for l in pl:
    x = xyz(l) - pc; cen.append(f'{l[:30]}{x[0]:8.3f}{x[1]:8.3f}{x[2]:8.3f}{l[54:]}')
open(f'{OUT}/pyrene_centered.pdb', 'w').write(''.join(cen) + 'END\n')

rng = random.Random(seed)
centres = []
tries = 0
while len(centres) < NP:
    tries += 1; assert tries < 100000, 'cannot place pyrene centres; lower the count or raise --rpyr'
    v = np.array([rng.uniform(-1, 1) for _ in range(3)])
    if np.linalg.norm(v) > 1: continue
    v = v * a.rpyr
    if all(np.linalg.norm(v - np.array(c)) >= a.pyr_sep for c in centres): centres.append(tuple(v))
chains = [f'{CH}/chain_s{(k % 4) + 1}_min.pdb' for k in range(NC)]
sizes = [len(pl)] * NP + [len(rec(c)) for c in chains]

inp = f'tolerance 2.0\nfiletype pdb\nseed {seed}\noutput {name}_raw.pdb\n'
for (x, y, z) in centres:
    ang = [rng.uniform(0, 2 * math.pi) for _ in range(3)]
    inp += f'\nstructure {OUT}/pyrene_centered.pdb\n  number 1\n  fixed {x:.3f} {y:.3f} {z:.3f} {ang[0]:.4f} {ang[1]:.4f} {ang[2]:.4f}\nend structure\n'
for c in chains:
    inp += f'\nstructure {c}\n  number 1\n  inside sphere 0. 0. 0. {R:.1f}\n'
    for (x, y, z) in centres: inp += f'  outside sphere {x:.3f} {y:.3f} {z:.3f} 12.\n'
    inp += 'end structure\n'
open(f'{OUT}/pack_{name}.inp', 'w').write(inp)
sh(f'packmol < pack_{name}.inp > pack_{name}.log 2>&1')
plog = open(f'{OUT}/pack_{name}.log').read()
check(os.path.exists(f'{OUT}/{name}_raw.pdb') and 'Success!' in plog, 'packmol finished successfully')
atoms = rec(f'{OUT}/{name}_raw.pdb')
check(len(atoms) == sum(sizes), f'raw pdb has {len(atoms)} atoms == {sum(sizes)}')
out = []; k = 0
for s in sizes: out += atoms[k:k + s] + ['TER\n']; k += s
open(f'{OUT}/{name}.pdb', 'w').write(''.join(out) + 'END\n')

# start-state check: every pyrene atom vs every chain atom
X = np.array([xyz(l) for l in atoms]); npy = 26 * NP
if NP and NC:
    from scipy.spatial.distance import cdist
    d = cdist(X[:npy], X[npy:]); print(f'  start: minimum pyrene-chain atom distance {d.min():.2f} A')
    check(d.min() > 6.0, 'no pyrene is in contact with a chain at the start (min atom distance > 6 A; Packmol rule is centre-to-atom >= 12 A)')
    pm = np.array([X[26 * i:26 * (i + 1)].mean(0) for i in range(NP)])
if NP > 1:
    from scipy.spatial.distance import cdist
    pm = np.array([X[26 * i:26 * (i + 1)].mean(0) for i in range(NP)])
    print('  start: minimum pyrene-pyrene centre distance %.1f A' % min(np.linalg.norm(pm[i] - pm[j]) for i in range(NP) for j in range(i)))
    pp = min(cdist(X[26 * i:26 * (i + 1)], X[26 * j:26 * (j + 1)]).min() for i in range(NP) for j in range(i))
    check(pp > 4.5, f'no pyrene-pyrene contact at the start (min atom distance {pp:.2f} A > 4.5 A)')
print(f'  start: packmol sphere R={R:.0f} A, {NP} pyrene, {NC} chains')

res_loads = ''.join(f'{r} = loadmol2 {P}/res_{r}.mol2\n' for r in ('MON', 'MNF', 'MNL', 'CPH', 'CPT'))
ff = 'source leaprc.gaff2\nsource leaprc.water.tip3p\n'
if NC: ff += f'loadamberparams {P}/monomer.frcmod\n{res_loads}'
if NP: ff += f'loadamberparams {PF}\nPYR = loadmol2 {PM}\n'
def solvate(buf):
    open(f'{OUT}/solv_{name}.in', 'w').write(ff + f'sys = loadpdb {name}.pdb\ncheck sys\ncharge sys\nsolvatebox sys TIP3PBOX {buf:.3f} iso\nsaveamberparm sys {name}.prmtop {name}.inpcrd\nquit\n')
    sh(f'tleap -f solv_{name}.in > leap_{name}.txt 2>&1')
    m_ = re.findall(r'Total bounding box for atom centers:\s+([\d.]+)', open(f'{OUT}/leap_{name}.txt').read())
    return float(m_[-1]) + 3.4 if m_ else None          # GROMACS box edge is about atom-centre box + 3.4 A (c8: 103.17 -> 106.5)
buf = a.buffer; edge = solvate(buf)
if a.box_edge and edge:
    for _ in range(2):
        buf = max(2.0, buf + (a.box_edge * 10 - edge) / 2); edge = solvate(buf)
    print(f'  box-edge target {a.box_edge:.2f} nm: buffer {buf:.2f} A -> tleap vdw box {edge / 10:.2f} nm')
log = open(f'{OUT}/leap_{name}.txt').read()
check('Errors = 0' in log, 'tleap errors = 0')
m = re.search(r'Total unperturbed charge:\s+(\S+)', log); check(m and abs(float(m.group(1))) < 1e-3, f'tleap net charge {m.group(1) if m else "?"}')
print(f'  tleap close-contact warnings: {len(re.findall("Close contact", log))}')
check(os.path.exists(f'{OUT}/{name}.prmtop'), 'prmtop written')

import parmed as pmd
s = pmd.load_file(f'{OUT}/{name}.prmtop', f'{OUT}/{name}.inpcrd')
s.save(f'{OUT}/{name}.top', overwrite=True); s.save(f'{OUT}/{name}.gro', overwrite=True)
top = open(f'{OUT}/{name}.top').read()
mol = re.search(r'\[ molecules \]\n(.*)', top, re.S).group(1)
print('  [ molecules ]:', ' | '.join(' '.join(l.split()) for l in mol.strip().split('\n') if not l.startswith(';')))
nwat = int(re.search(r'WAT\s+(\d+)', mol).group(1))
box = open(f'{OUT}/{name}.gro').read().split('\n')[-2].split()
vol = np.prod([float(b) for b in box[:3]])
print(f'  box (nm): {box}  atoms: {len(s.atoms)}  waters: {nwat}')
if NP: print(f'  pyrene concentration {NP / (vol * 0.6022) * 1000:.1f} mM; chain concentration {NC / (vol * 0.6022) * 1000:.1f} mM')
check(len(s.atoms) == 26 * NP + 556 * NC + 3 * nwat, f'atom count {len(s.atoms)} == solute {26 * NP + 556 * NC} + 3 x {nwat} waters')
if NP:
    q = np.array([at.charge for at in s.atoms[:26 * NP]]).reshape(NP, 26)
    check(np.abs(q.sum(1)).max() < 1e-4 and abs(q[:, 16:].max() - 0.17) < 0.01, f'pyrene charges corrected (max H {q[:, 16:].max():.3f}, max |sum| {np.abs(q.sum(1)).max():.1e})')
print('RESULT:', 'ALL CHECKS PASSED' if ok else 'FAILED'); sys.exit(0 if ok else 1)
