#!/usr/bin/env python3
"""Per-frame geometry for the dashboard (scripts/build_site.py), from a whole-molecule FULL-SYSTEM xtc (gmx trjconv -pbc whole -dt 500).

  site_metrics.py --gro md_X.gro --xtc whole500.xtc --npyr NP --nchains NC --out metrics.csv

Atom order (build_box.py / build_many.py): NP pyrene of 26 atoms, then NC chains of 556 atoms, then TIP3P water (resname WAT, oxygen 'O').
Per frame (heavy atoms, minimum image; same definitions as render_*_movie.py and analyze_many.py):
  waters_total        waters whose O is within --shell nm (0.45) of any pyrene heavy atom (union over pyrene)
  waters_pyrK         same, for pyrene K alone
  rg_chainK           chain K radius of gyration (nm; masses guessed from atom names, matches gmx gyrate to ~0.002 nm)
  chains_on_pyrK      chains with a heavy atom within --contact nm (0.4) of pyrene K
  mindist_pyrK        min pyrene K - any chain heavy-atom distance (nm)
  pyr_clusters / largest_pyr_cluster     components of the pyrene-pyrene contact graph
  chain_clusters / largest_chain_cluster components of the chain-chain contact graph
  pyr_free            pyrene with no chain contact
  mindist_chainK      min chain K - any pyrene heavy-atom distance (nm)
"""
import argparse, warnings
import numpy as np
import MDAnalysis as mda
from MDAnalysis.lib.distances import capped_distance, distance_array
warnings.filterwarnings('ignore')

ap = argparse.ArgumentParser()
ap.add_argument('--gro', required=True); ap.add_argument('--xtc', required=True); ap.add_argument('--out', required=True)
ap.add_argument('--npyr', type=int, required=True); ap.add_argument('--nchains', type=int, required=True)
ap.add_argument('--sasa', action='store_true', help='also pyrene SASA coverage (Shrake-Rupley on pyrene atoms only; 0.14 nm probe, 960 dots, Bondi radii)')
ap.add_argument('--shell', type=float, default=0.45); ap.add_argument('--contact', type=float, default=0.4)
a = ap.parse_args()
NP, NC, P, C = a.npyr, a.nchains, 26, 556
u = mda.Universe(a.gro, a.xtc)
names = u.atoms.names; resn = u.atoms.resnames
heavy = np.array([not n.startswith('H') for n in names])
mol = np.full(len(u.atoms), -1)
for i in range(NP): mol[i * P:(i + 1) * P] = i
for k in range(NC): mol[NP * P + k * C: NP * P + (k + 1) * C] = NP + k
sol_h = np.where(heavy & (mol >= 0))[0]; sol_mol = mol[sol_h]
pyr_h = sol_h[sol_mol < NP]; pyr_of = sol_mol[sol_mol < NP]
wat_o = np.where((resn == 'WAT') & (names == 'O'))[0]
pyr_hi = [pyr_h[pyr_of == i] for i in range(NP)]
ch_h = [sol_h[sol_mol == NP + k] for k in range(NC)]
chains = [u.atoms[NP * P + k * C: NP * P + (k + 1) * C] for k in range(NC)]
nm = NP + NC

def components(adj):
    n = len(adj); lab = -np.ones(n, int); c = 0
    for s in range(n):
        if lab[s] >= 0: continue
        st = [s]; lab[s] = c
        while st:
            v = st.pop()
            for w in np.where(adj[v] & (lab < 0))[0]: lab[w] = c; st.append(w)
        c += 1
    return lab

# ---- pyrene-only Shrake-Rupley: same definition as analyze_many.py (1 - SASA(pyr | occluders)/SASA(pyr alone)), but dots only on pyrene atoms
RAD = {'C': 1.70, 'H': 1.20, 'N': 1.55, 'O': 1.52, 'S': 1.80}; PROBE = 1.4                       # Angstrom (Bondi radii, as gmx sasa)
rad = np.array([RAD.get(n[0], 1.70) for n in names])
gi = np.arange(960) + 0.5; phi = np.arccos(1 - 2 * gi / 960); th = np.pi * (1 + 5 ** 0.5) * gi
DOTS = np.stack([np.cos(th) * np.sin(phi), np.sin(th) * np.sin(phi), np.cos(phi)], 1)
sol_all = np.where(mol >= 0)[0]; sol_all_mol = mol[sol_all]
pyr_all = [np.where(mol == i)[0] for i in range(NP)]
def pyr_sasa(X, box):
    # returns per pyrene: SASA alone, | chains + self, | all solute   (A^2)
    out = np.zeros((NP, 3)); L = box[:3]
    for i in range(NP):
        pa = pyr_all[i]
        pr = capped_distance(X[pa], X[sol_all], 2 * 1.80 + 2 * PROBE, box=box, return_distances=False)
        for ai, a_ in enumerate(pa):
            J = sol_all[pr[pr[:, 0] == ai, 1]]; J = J[J != a_]
            R = rad[a_] + PROBE; P = X[a_] + R * DOTS
            d = X[J][None, :, :] - P[:, None, :]; d -= L * np.round(d / L)
            hit = (d * d).sum(2) < (rad[J] + PROBE)[None, :] ** 2
            cat = sol_all_mol[np.searchsorted(sol_all, J)]
            own, chn = cat == i, cat >= NP
            b_alone = hit[:, own].any(1); b_poly = b_alone | hit[:, chn].any(1); b_tot = hit.any(1)
            area = 4 * np.pi * R * R / len(DOTS)
            out[i] += area * np.array([(~b_alone).sum(), (~b_poly).sum(), (~b_tot).sum()])
    return out

cols = ['time_ns', 'waters_total'] + [f'waters_pyr{i+1}' for i in range(NP)] + [f'rg_chain{k+1}' for k in range(NC)] \
     + [f'chains_on_pyr{i+1}' for i in range(NP)] + [f'mindist_pyr{i+1}' for i in range(NP)] \
     + ['pyr_clusters', 'largest_pyr_cluster', 'chain_clusters', 'largest_chain_cluster', 'pyr_free'] \
     + ([f'mindist_chain{k+1}' for k in range(NC)] if (NP and NC) else []) \
     + ([f'covpoly_pyr{i+1}' for i in range(NP)] * bool(NC) + [f'covtot_pyr{i+1}' for i in range(NP)] + (['covpoly_mean'] if NC else []) + ['covtot_mean'] if (a.sasa and NP) else [])
rows = []
for ts in u.trajectory:
    X = u.atoms.positions; box = ts.dimensions; r = [ts.time / 1000.0]
    if NP and len(wat_o):
        pr = capped_distance(X[pyr_h], X[wat_o], a.shell * 10, box=box, return_distances=False)
        r.append(len(np.unique(pr[:, 1])))
        r += [len(np.unique(pr[pyr_of[pr[:, 0]] == i, 1])) for i in range(NP)]
    else: r += [np.nan] + [np.nan] * NP
    r += [c.radius_of_gyration() / 10.0 for c in chains]
    # molecule contact graph from one capped self-distance search over solute heavy atoms
    pairs = capped_distance(X[sol_h], X[sol_h], a.contact * 10, box=box, return_distances=False)
    mi, mj = sol_mol[pairs[:, 0]], sol_mol[pairs[:, 1]]; keep = mi != mj
    A = np.zeros((nm, nm), bool); A[mi[keep], mj[keep]] = True; A |= A.T
    on = A[:NP, NP:].sum(1) if NC else np.zeros(NP)
    dpc = np.array([[distance_array(X[pyr_hi[i]], X[ch_h[k]], box=box).min() / 10.0 for k in range(NC)] for i in range(NP)]).reshape(NP, NC)
    r += list(on) + (list(dpc.min(1)) if NC else [np.nan] * NP)
    if NP:
        lp = components(A[:NP, :NP]); r += [lp.max() + 1, np.bincount(lp).max()]
    else: r += [np.nan, np.nan]
    if NC:
        lc = components(A[NP:, NP:]); r += [lc.max() + 1, np.bincount(lc).max()]
    else: r += [np.nan, np.nan]
    r.append(int((on == 0).sum()) if (NP and NC) else np.nan)
    r += list(dpc.min(0)) if (NP and NC) else []
    if a.sasa and NP:
        sa = pyr_sasa(X, box); cp = 100 * (1 - sa[:, 1] / sa[:, 0]); ct = 100 * (1 - sa[:, 2] / sa[:, 0])
        r += (list(cp) if NC else []) + list(ct) + ([cp.mean()] if NC else []) + [ct.mean()]
    rows.append(r)
with open(a.out, 'w') as f:
    f.write(','.join(cols) + '\n')
    for r in rows: f.write(','.join('nan' if (isinstance(v, float) and np.isnan(v)) else (f'{v:.4f}' if isinstance(v, (float, np.floating)) else str(int(v))) for v in r) + '\n')
print(f'{a.out}: {len(rows)} frames, {rows[0][0]:.1f}-{rows[-1][0]:.1f} ns')
