#!/usr/bin/env python3
"""Build a ROMP chain with a chosen cis/trans pattern at every backbone C=C junction, from tleap residues.

  build_chain.py --params params/monomer_c1 --out chains/c1_dp10 --tag s1 --dp 10 --pattern CCCCTTTTC
  build_chain.py ... --seed 7                      random 50/50 per junction (reproducible)
  build_chain.py ... --like-old chains/dp10/chain_s1_min.pdb   copy the cis/trans pattern of an older chain

Pattern letters run from the CPH end: junction k joins monomer k and k+1 (dp-1 junctions). C = cis (Z, dihedral ~0),
T = trans (E, ~180). Dihedral = ring carbon - tail C = head C - ring carbon of the next monomer.

Method (same as the first chain set, generalised to any DP): tleap `sequence {CPH MNF MON... MNL CPT}`; each junction is set by
rotating everything downstream rigidly about the C=C axis; tleap prmtop; vacuum sander minimisation; then validation.
Needs tleap, sander, ambpdb (AmberTools) and numpy. Atom order is identical for every pattern, so one topology serves all patterns.
Exit code 1 if any validation check fails.
"""
import argparse, os, random, re, subprocess, sys
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument('--params', default='params/monomer_c1'); ap.add_argument('--out', required=True); ap.add_argument('--tag', required=True)
ap.add_argument('--dp', type=int, default=10)
g = ap.add_mutually_exclusive_group(required=True)
g.add_argument('--pattern'); g.add_argument('--seed', type=int); g.add_argument('--like-old')
ap.add_argument('--old-atoms', default='C7,C8,C4,C3', help='names (ring C, tail C, head C, ring C) used to read an old chain')
ap.add_argument('--frcmod', default=None); ap.add_argument('--maxcyc', type=int, default=10000)
a = ap.parse_args()
assert a.dp >= 3, 'dp must be >= 3 (CPH MNF ... MNL CPT)'
P = os.path.abspath(a.params); OUT = os.path.abspath(a.out); os.makedirs(OUT, exist_ok=True)
frcmod = a.frcmod or os.path.join(P, 'monomer.frcmod')
run = lambda cmd: subprocess.run(cmd, shell=True, cwd=OUT, capture_output=True, text=True)

def read_mol2(path):
    raw = open(path).read(); sec = {p.split('\n', 1)[0].strip(): p.split('\n', 1)[1].strip().split('\n') for p in re.split(r'@<TRIPOS>', raw)[1:]}
    at = [l.split() for l in sec['ATOM']]; bo = [(int(l.split()[1]) - 1, int(l.split()[2]) - 1) for l in sec['BOND']]
    return [x[1] for x in at], np.array([[float(v) for v in x[2:5]] for x in at]), bo, [l.split()[0] for l in sec['HEADTAIL']]
names, ref_xyz, bonds, ht = read_mol2(os.path.join(P, 'res_MON.mol2'))
nb = {i: [] for i in range(len(names))}
for i, j in bonds: nb[i].append(j); nb[j].append(i)
head, tail = ht
hi, ti = names.index(head), names.index(tail)
heavy_nb = lambda c: [n for n in nb[c] if not names[n].startswith('H')]
assert len(heavy_nb(hi)) == 1 and len(heavy_nb(ti)) == 1, 'head/tail carbons should each have one heavy neighbour inside the residue'
hring, tring = names[heavy_nb(hi)[0]], names[heavy_nb(ti)[0]]

njunc = a.dp - 1
def dihedral(p):
    p = [np.asarray(q, float) for q in p]
    b1, b2, b3 = p[1] - p[0], p[2] - p[1], p[3] - p[2]
    n1, n2 = np.cross(b1, b2), np.cross(b2, b3); m = np.cross(n1, b2 / np.linalg.norm(b2))
    return float(np.degrees(np.arctan2(np.dot(m, n2), np.dot(n1, n2))))
def read_pdb(path):
    L = open(path).read().splitlines()
    atoms = [dict(res=int(l[22:26]), name=l[12:16].strip(), xyz=np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])) for l in L if l.startswith(('ATOM', 'HETATM'))]
    for k, x in enumerate(atoms): x['i'] = k
    return L, atoms

if a.pattern: pat = a.pattern.upper()
elif a.seed is not None:
    rr = random.Random(a.seed); pat = ''.join(rr.choice('CT') for _ in range(njunc))
else:
    _, oa = read_pdb(a.like_old); n1, n2, n3, n4 = a.old_atoms.split(',')
    d = {(x['name'], x['res']): x['xyz'] for x in oa}; pat = ''
    for r in range(2, 2 + njunc):
        pat += 'C' if abs(dihedral([d[(n1, r)], d[(n2, r)], d[(n3, r + 1)], d[(n4, r + 1)]])) < 90 else 'T'
assert len(pat) == njunc and set(pat) <= set('CT'), f'pattern must be {njunc} letters of C/T, got {pat!r}'
print(f'pattern ({njunc} junctions): {pat}   cis {pat.count("C")}  trans {pat.count("T")}')

seq = 'CPH MNF ' + ' '.join(['MON'] * (a.dp - 2)) + ' MNL CPT'
loads = ''.join(f'{r} = loadmol2 {os.path.join(P, "res_" + r + ".mol2")}\n' for r in ('MON', 'MNF', 'MNL', 'CPH', 'CPT'))
def leap(inp, tag):
    open(os.path.join(OUT, tag + '.in'), 'w').write(inp); run(f'tleap -f {tag}.in > {tag}.log 2>&1')
    log = open(os.path.join(OUT, tag + '.log')).read(); assert 'Errors = 0' in log, f'tleap errors, see {tag}.log'
raw = f'chain_{a.tag}_raw.pdb'
leap(f'source leaprc.gaff2\nloadamberparams {frcmod}\n{loads}chain = sequence {{ {seq} }}\ncheck chain\ncharge chain\nsavepdb chain {raw}\nquit\n', 'build_' + a.tag)
L, atoms = read_pdb(os.path.join(OUT, raw)); nat = len(atoms); print('atoms', nat)

byres = {}
for x in atoms: byres.setdefault(x['res'], {})[x['name']] = x['i']
cur = {x['i']: x['xyz'].copy() for x in atoms}
def rot(points, p0, axis, th):
    k = axis / np.linalg.norm(axis); v = points - p0; c, s = np.cos(th), np.sin(th)
    return p0 + v * c + np.cross(k, v) * s + np.outer(v @ k, k) * (1 - c)
def junction_dih(r):
    return dihedral([cur[byres[r][tring]], cur[byres[r][tail]], cur[byres[r + 1][head]], cur[byres[r + 1][hring]]])
for k, ch in enumerate(pat):
    r = 2 + k; target = 0.0 if ch == 'C' else 180.0
    delta = np.radians(target - junction_dih(r))
    p0 = cur[byres[r][tail]]; axis = cur[byres[r + 1][head]] - p0
    idx = [x['i'] for x in atoms if x['res'] > r]; save = {i: cur[i].copy() for i in idx}
    for sgn in (1, -1):
        for i, q in zip(idx, rot(np.array([save[i] for i in idx]), p0, axis, sgn * delta)): cur[i] = q
        d1 = junction_dih(r)
        if min(abs(d1 - target), 360 - abs(d1 - target)) < 1.0: break
    else: sys.exit(f'could not set junction {k + 1}')
ai = 0; new_lines = []
for l in L:
    if l.startswith(('ATOM', 'HETATM')):
        p = cur[ai]; l = l[:30] + f'{p[0]:8.3f}{p[1]:8.3f}{p[2]:8.3f}' + l[54:]; ai += 1
    new_lines.append(l)
# tleap's savepdb writes a TER after every residue; loadpdb would treat each residue as its own chain and make no
# inter-residue bonds. Keep one TER at the very end only.
new_lines = [l for l in new_lines if not l.startswith('TER')]
new_lines.insert(max(i for i, l in enumerate(new_lines) if l.startswith(('ATOM', 'HETATM'))) + 1, 'TER')
flip_pdb = f'chain_{a.tag}.pdb'; open(os.path.join(OUT, flip_pdb), 'w').write('\n'.join(new_lines) + '\n')

leap(f'source leaprc.gaff2\nloadamberparams {frcmod}\n{loads}chain = loadpdb {flip_pdb}\ncheck chain\ncharge chain\nsaveamberparm chain chain_{a.tag}.prmtop chain_{a.tag}.inpcrd\nquit\n', 'flip_' + a.tag)
open(os.path.join(OUT, 'min_chain.in'), 'w').write(f'vacuum minimization\n&cntrl\n imin=1, maxcyc={a.maxcyc}, ncyc={a.maxcyc // 5},\n ntb=0, igb=0, cut=999.,\n/\n')
r = run(f'sander -O -i min_chain.in -o min_{a.tag}.out -p chain_{a.tag}.prmtop -c chain_{a.tag}.inpcrd -r chain_{a.tag}_min.rst')
assert os.path.exists(os.path.join(OUT, f'chain_{a.tag}_min.rst')), 'sander failed: ' + r.stderr[-300:]
run(f'ambpdb -p chain_{a.tag}.prmtop -c chain_{a.tag}_min.rst > chain_{a.tag}_min.pdb')
mt = open(os.path.join(OUT, f'min_{a.tag}.out')).read()
m = re.search(r'FINAL RESULTS.*?NSTEP\s+ENERGY\s+RMS\s+GMAX.*?\n\s+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)', mt, re.S)
print('sander final: nstep %s  E %s  RMS %s  Gmax %s' % m.groups() if m else 'sander final: see out file')

ok = True
def check(cond, msg):
    global ok
    print(('PASS ' if cond else 'FAIL ') + msg); ok &= bool(cond)
L2, A2 = read_pdb(os.path.join(OUT, f'chain_{a.tag}_min.pdb'))
check(len(A2) == nat, f'atom count {len(A2)} == {nat}')
check([(x['name'], x['res']) for x in A2] == [(x['name'], x['res']) for x in atoms], 'atom names and residue numbers identical to the tleap build (one topology serves every pattern)')
br = {}
for x in A2: br.setdefault(x['res'], {})[x['name']] = x['xyz']
got = ''; worst = 0
for k in range(njunc):
    r = 2 + k; d = dihedral([br[r][tring], br[r][tail], br[r + 1][head], br[r + 1][hring]])
    got += 'C' if abs(d) < 90 else 'T'; worst = max(worst, min(abs(d), 180 - abs(d)))
check(got == pat, f'pattern after minimisation {got} == requested {pat}')
nbond = [float(np.linalg.norm(br[r][tail] - br[r + 1][head])) for r in range(2, a.dp + 1)]
check(max(nbond) < 1.5, f'inter-residue C-C bonds present ({min(nbond):.2f}-{max(nbond):.2f} A)')
check(worst < 30, f'worst junction deviation from planar {worst:.1f} deg (< 30)')
cc = [float(np.linalg.norm(br[2 + k][tail] - br[3 + k][head])) for k in range(njunc)]
check(all(1.30 < c < 1.40 for c in cc), f'junction C=C lengths {min(cc):.3f}-{max(cc):.3f} A (expect 1.30-1.40)')
bad = 0
for r in range(1, a.dp + 3):
    rn = br[r]
    for i, j in bonds:
        if names[i] in rn and names[j] in rn and not (0.9 < np.linalg.norm(rn[names[i]] - rn[names[j]]) < 1.8): bad += 1
check(bad == 0, f'{bad} template bonds outside 0.9-1.8 A')
sgn = lambda xyz, c, ns: np.sign(np.linalg.det(np.array([xyz[n] - xyz[c] for n in ns])))
cent = [(c, nb[c][:3]) for c in range(len(names)) if len(nb[c]) == 4]
refs = {c: sgn(ref_xyz, c, ns) for c, ns in cent}; flipped = tot = 0
for r in range(3, a.dp + 1):
    xyz = np.array([br[r][n] for n in names])
    for c, ns in cent: tot += 1; flipped += sgn(xyz, c, ns) != refs[c]
check(flipped == 0, f'stereo retained at {tot} sp3 centres in internal residues ({flipped} inverted)')
Xall = np.array([x['xyz'] for x in A2]); D2 = np.linalg.norm(Xall[:, None] - Xall[None], axis=2); np.fill_diagonal(D2, 9)
check(D2.min() > 0.7, f'no atom pair closer than 0.7 A (min {D2.min():.2f})')
hv = [x for x in A2 if not x['name'].startswith('H')]; X = np.array([x['xyz'] for x in hv]); R = np.array([x['res'] for x in hv])
mind = np.linalg.norm(X[:, None] - X[None], axis=2)[np.abs(R[:, None] - R[None]) >= 2].min()
check(mind > 2.0, f'closest heavy-atom contact between non-adjacent residues {mind:.2f} A (need > 2.0)')
print('extent (A):', np.round(np.ptp(Xall, axis=0), 1), ' charge:', re.findall(r'Total unperturbed charge:\s+(\S+)', open(os.path.join(OUT, 'build_' + a.tag + '.log')).read()))
print('RESULT:', 'ALL CHECKS PASSED' if ok else 'FAILED'); sys.exit(0 if ok else 1)
