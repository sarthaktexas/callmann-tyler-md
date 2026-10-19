"""Build the unbound starting solute for umbrella pattern N (chain_sN_min.pdb + pyrene).

Writes <outdir>/solute_unboxed.gro (pyrene = atoms 1-26, chain = 27-582, same order as boxes/sys_1to1.top).
Pyrene is placed 3.1 nm (COM-COM) from the CORE group COM, in the direction (of 2000 random ones) that keeps it
farthest from every chain atom, and must be >= 1.2 nm from all chain atoms, including periodic images (7 nm box). The same recipe is used for all four patterns.
Usage: python build_start.py <pattern 1-4> <outdir> [box_nm=7.0]
"""
import sys, os, numpy as np
N = int(sys.argv[1]); out = sys.argv[2]; box = float(sys.argv[3]) if len(sys.argv) > 3 else 7.0
HERE = os.path.dirname(os.path.abspath(__file__)); REPO = os.path.abspath(os.path.join(HERE, '..', '..'))
MASS = {'H': 1.008, 'C': 12.011, 'N': 14.007, 'O': 15.999}
def read_pdb(fn):
    rec = []
    for l in open(fn):
        if l.startswith(('ATOM', 'HETATM')):
            rec.append((l[12:16].strip(), l[17:20].strip(), int(l[22:26]), np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])]) / 10.0))
    return rec
pyr = read_pdb(os.path.join(REPO, 'params/pyrene/pyrene_gaff2.pdb')); chain = read_pdb(os.path.join(REPO, f'chains/dp10/chain_s{N}_min.pdb'))
assert len(pyr) == 26 and len(chain) == 556, (len(pyr), len(chain))
P = np.array([a[3] for a in pyr]); C = np.array([a[3] for a in chain])
core = [int(x) for x in open(os.path.join(HERE, 'core_heavy_atoms.txt')).read().split()]
assert len(core) == 120 and min(core) > 26 and max(core) <= 582
core_xyz = C[[i - 27 for i in core]]; core_m = np.array([MASS[chain[i - 27][0][0]] for i in core])
core_com = (core_xyz * core_m[:, None]).sum(0) / core_m.sum()
pm = np.array([MASS[a[0][0]] for a in pyr]); P0 = P - (P * pm[:, None]).sum(0) / pm.sum()   # pyrene centred on its COM
rng = np.random.default_rng(N); best = None
for _ in range(2000):
    d = rng.normal(size=3); d /= np.linalg.norm(d)
    com = core_com + 3.1 * d; X = P0 + com
    diff = X[:, None, :] - C[None, :, :]; diff -= box * np.round(diff / box)      # minimum image
    dmin = np.sqrt((diff ** 2).sum(2)).min()
    if best is None or dmin > best[0]: best = (dmin, X, d)
dmin, X, d = best
print('pattern %d: min pyrene-chain atom distance %.2f nm, direction %s' % (N, dmin, np.round(d, 3)))
assert dmin >= 1.2, 'pyrene too close to chain; lengthen the 3.1 nm start distance or use a larger box'
allxyz = np.vstack([X, C]); shift = box / 2 - allxyz.mean(0); allxyz += shift
atoms = pyr + chain
with open(os.path.join(out, 'solute_unboxed.gro'), 'w') as f:
    f.write('pyrene + chain_s%d, unbound start\n%d\n' % (N, len(atoms)))
    for i, (a, xyz) in enumerate(zip(atoms, allxyz)):
        f.write('%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n' % (a[2] % 100000, a[1], a[0], (i + 1) % 100000, *xyz))
    f.write('%10.5f%10.5f%10.5f\n' % (box, box, box))
# Rg-like sanity: extent of the chain (it starts extended)
print('chain extent (max pairwise-ish): %.2f nm, box %.1f nm' % (np.ptp(allxyz[26:], axis=0).max(), box))
