#!/usr/bin/env python3
"""Load pyrene into a PRE-ASSEMBLED (chain-only, equilibrated) polymer box: the 'pre-formed micelle' experiment.

  python3 scripts/load_pyrene.py --gro boxes_many/md_many_p0_c16.gro --n-chain 16 --n-pyr 8 \
        --template boxes_many/many_p8_c16.gro --top boxes_many/many_p8_c16.top --out boxes_many/many_p8_c16_pre

Takes the final frame of a chain-only run (n_chain DP=10 chains + TIP3P), places n_pyr pyrene in the WATER with random
orientation, removes overlapping waters and writes <out>.gro / <out>.top with the same atom order as the co-assembled
box (all pyrene first, then chains, then water), so the co-assembled topology (moleculetypes PYR, system1, WAT) is reused
with only the WAT count changed and analyze_many.py works unchanged.

Placement rules (minimum image, orthorhombic box):
  * every pyrene atom >= --min-chain nm (default 1.5) from every chain atom  -> unbound start, by construction
  * pyrene centres >= --min-pyr nm apart (default 3.0, as build_many.py)
  * a water molecule is removed if any of its atoms is < --water-cut nm (default 0.30) from any pyrene atom
Velocities are dropped (nvt.mdp has gen-vel = yes, so the loaded box starts from fresh velocities).
Standard library only, so it runs with any python3 on ARC.
"""
import argparse, math, random, re, sys

ap = argparse.ArgumentParser()
ap.add_argument('--gro', required=True, help='final frame of the chain-only run (e.g. md_many_p0_c16.gro)')
ap.add_argument('--n-chain', type=int, required=True); ap.add_argument('--n-pyr', type=int, default=8)
ap.add_argument('--template', required=True, help='gro whose first 26 atoms are one PYR (e.g. many_p8_c16.gro)')
ap.add_argument('--top', required=True, help='co-assembled topology with PYR, system1, WAT moleculetypes')
ap.add_argument('--out', required=True, help='output prefix (writes <out>.gro, <out>.top)')
ap.add_argument('--min-chain', type=float, default=1.5); ap.add_argument('--min-pyr', type=float, default=3.0)
ap.add_argument('--water-cut', type=float, default=0.30); ap.add_argument('--seed', type=int, default=4816)
a = ap.parse_args()
NAT_CH, NAT_PY = 556, 26
ok = True
def check(c, m):
    global ok; print(('PASS ' if c else 'FAIL ') + m); ok = ok and bool(c)

def read_gro(p):
    L = open(p).read().splitlines()
    n = int(L[1]); at = []
    for l in L[2:2 + n]:
        at.append((int(l[0:5]), l[5:10].strip(), l[10:15].strip(), float(l[20:28]), float(l[28:36]), float(l[36:44])))
    box = [float(x) for x in L[2 + n].split()]
    return at, box

atoms, box = read_gro(a.gro)
assert len(box) == 3 or all(abs(x) < 1e-6 for x in box[3:]), 'triclinic box not supported'
B = box[:3]
nch = NAT_CH * a.n_chain
chain, water = atoms[:nch], atoms[nch:]
check(len(water) % 3 == 0 and all(w[1] == 'WAT' for w in water), f'input: {a.n_chain} chains ({nch} atoms) then {len(water) // 3} WAT')
check(not any(c[1] == 'WAT' for c in chain), 'no WAT inside the chain block')
tmpl, _ = read_gro(a.template)
pyr = tmpl[:NAT_PY]
check(all(p[1] == 'PYR' for p in pyr), 'template: first 26 atoms are PYR')

# same atom names/order for the chains as in the co-assembled box (so its system1 moleculetype applies)
tall = tmpl
ref_chain = [t[2] for t in tall[NAT_PY * a.n_pyr:NAT_PY * a.n_pyr + NAT_CH]]
check(all([c[2] for c in chain[k * NAT_CH:(k + 1) * NAT_CH]] == ref_chain for k in range(a.n_chain)), 'chain atom names match the template system1 order')

def d2(p, q):
    s = 0.0
    for i in range(3):
        x = p[i] - q[i]; x -= B[i] * round(x / B[i]); s += x * x
    return s

C = [(c[3], c[4], c[5]) for c in chain]
pc = [sum(p[3 + i] for p in pyr) / NAT_PY for i in range(3)]
P0 = [(p[3] - pc[0], p[4] - pc[1], p[5] - pc[2]) for p in pyr]
rmax = max(math.sqrt(sum(x * x for x in p)) for p in P0)

def rot(rng):
    u1, u2, u3 = rng.random(), rng.random(), rng.random()          # uniform random quaternion (Shoemake)
    q = (math.sqrt(1 - u1) * math.sin(2 * math.pi * u2), math.sqrt(1 - u1) * math.cos(2 * math.pi * u2),
         math.sqrt(u1) * math.sin(2 * math.pi * u3), math.sqrt(u1) * math.cos(2 * math.pi * u3))
    x, y, z, w = q
    return ((1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)),
            (2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)),
            (2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)))

rng = random.Random(a.seed)
placed, cents, tries = [], [], 0
mc2, mp2 = a.min_chain ** 2, a.min_pyr ** 2
while len(placed) < a.n_pyr:
    tries += 1
    if tries > 20000:
        sys.exit(f'FAIL could only place {len(placed)} pyrene; lower --min-chain/--min-pyr')
    c = [rng.random() * B[i] for i in range(3)]
    if any(d2(c, q) < mp2 for q in cents): continue
    R = rot(rng)
    X = [tuple(c[i] + sum(R[i][j] * p[j] for j in range(3)) for i in range(3)) for p in P0]
    near = [x for x in C if d2(c, x) < (a.min_chain + rmax) ** 2]   # only these chain atoms can violate the rule
    if any(d2(xp, x) < mc2 for xp in X for x in near): continue
    placed.append(X); cents.append(c)
print(f'  placed {a.n_pyr} pyrene after {tries} trial centres (box {B[0]:.3f} x {B[1]:.3f} x {B[2]:.3f} nm)')

PY = [x for X in placed for x in X]
dmin_ch = math.sqrt(min(d2(p, x) for p in PY for x in C))
dmin_pp = math.sqrt(min(d2(cents[i], cents[j]) for i in range(len(cents)) for j in range(i)))
check(dmin_ch >= a.min_chain - 1e-9, f'minimum pyrene-chain atom distance {dmin_ch:.2f} nm (>= {a.min_chain})')
check(dmin_pp >= a.min_pyr - 1e-9, f'minimum pyrene-pyrene centre distance {dmin_pp:.2f} nm (>= {a.min_pyr})')

wc2, reach2 = a.water_cut ** 2, (a.water_cut + rmax + 0.2) ** 2
keep = []
for k in range(0, len(water), 3):
    W = water[k:k + 3]; o = (W[0][3], W[0][4], W[0][5])
    near = [i for i, cc in enumerate(cents) if d2(o, cc) < reach2]
    if near and any(d2((w[3], w[4], w[5]), p) < wc2 for w in W for i in near for p in placed[i]):
        continue
    keep.append(W)
nw0, nw = len(water) // 3, len(keep)
print(f'  waters: {nw0} -> {nw} ({nw0 - nw} removed, {(nw0 - nw) / a.n_pyr:.1f} per pyrene)')
check(0 < nw0 - nw < 40 * a.n_pyr, 'number of removed waters is plausible for pyrene-sized cavities')

wrap = lambda x: tuple(x[i] - B[i] * math.floor(x[i] / B[i]) for i in range(3))
lines, n = [], 0
def emit(resnr, resname, name, x):
    global n
    n += 1; lines.append(f'{resnr % 100000:5d}{resname:<5s}{name:>5s}{n % 100000:5d}{x[0]:8.3f}{x[1]:8.3f}{x[2]:8.3f}')
for k, X in enumerate(placed):
    cw = wrap(cents[k]); sh = [cw[i] - cents[k][i] for i in range(3)]   # keep each pyrene whole, centre inside the box
    for p, x in zip(pyr, X): emit(k + 1, 'PYR', p[2], (x[0] + sh[0], x[1] + sh[1], x[2] + sh[2]))
for c in chain: emit(c[0] + a.n_pyr, c[1], c[2], (c[3], c[4], c[5]))
r0 = chain[-1][0] + a.n_pyr if chain else a.n_pyr
for j, W in enumerate(keep):
    for w in W: emit(r0 + j + 1, 'WAT', w[2], (w[3], w[4], w[5]))
title = f'pre-assembled {a.n_chain} chains from {a.gro} + {a.n_pyr} pyrene loaded (load_pyrene.py seed {a.seed})'
open(a.out + '.gro', 'w').write(title + '\n' + f'{n:5d}\n' + '\n'.join(lines) + '\n' + ''.join(f'{b:10.5f}' for b in B) + '\n')

top = open(a.top).read()
mol = re.search(r'\[ molecules \]\n(?:;[^\n]*\n)*', top)
body = top[mol.end():]
check(re.search(rf'^PYR\s+{a.n_pyr}\s*$', body, re.M) and re.search(rf'^system1\s+{a.n_chain}\s*$', body, re.M),
      f'template topology has PYR {a.n_pyr} and system1 {a.n_chain}')
newbody = re.sub(r'^(WAT\s+)\d+', lambda m: f'{m.group(1)}{nw}', body, flags=re.M)
check(re.search(rf'^WAT\s+{nw}\s*$', newbody, re.M), f'output topology WAT count set to {nw}')
open(a.out + '.top', 'w').write(top[:mol.end()] + newbody)
check(n == NAT_PY * a.n_pyr + nch + 3 * nw, f'atom count {n} == {NAT_PY * a.n_pyr} pyrene + {nch} chain + 3 x {nw} water')
vol = B[0] * B[1] * B[2]
print(f'  pyrene concentration {a.n_pyr / (vol * 0.6022) * 1000:.1f} mM; chain concentration {a.n_chain / (vol * 0.6022) * 1000:.1f} mM')
print('RESULT:', 'ALL CHECKS PASSED' if ok else 'FAILED'); sys.exit(0 if ok else 1)
