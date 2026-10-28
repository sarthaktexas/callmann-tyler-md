#!/usr/bin/env python3
"""Mechanism analyses for one box, from the dashboard's whole-molecule 500 ps trajectory (site/work/<name>/whole500.xtc).

  mechanism_analysis.py --name many_p8_c16 --gro boxes_many/md_many_p8_c16.gro --npyr 8 --nchains 16 [--xtc ...] [--out analysis/mechanism]

Writes <out>/<name>_ts.csv (per frame), <name>_profiles.csv (radial densities, t >= 20 ns), <name>_pairs.npz (pyrene x chain contact
matrix per frame) and <name>_contacts.csv (contacting chain atoms by group, t >= 20 ns).
Definitions (heavy atoms, minimum image unless stated):
  aggregate     connected component of the chain-chain contact graph (min heavy-atom distance < 0.4 nm); molecules unwrapped around it
  agg_size / agg_rg / agg_kappa2   chains in the largest aggregate, its heavy-atom radius of gyration (nm), relative shape anisotropy
                (0 sphere ... 1 rod)
  depth_pyrK    distance of pyrene K's centroid from the centre of the aggregate it touches, divided by that aggregate's
                R = sqrt(5/3) Rg (0 = centre, ~1 = surface); NaN if the pyrene touches no chain
  stack_imide_pyrK  imide rings (C8 C9 N1 C20 C21) stacked on pyrene K: ring planes within 30 deg, imide centroid 3.0-4.2 A from the pyrene
                plane and within 4.0 A of the pyrene centroid in-plane
  pp_stacked    pyrene pairs stacked on each other (same geometric test, pyrene-pyrene)
  hb_inter / hb_intra  glucose O-H...O hydrogen bonds between different chains / within one chain (O...O < 3.5 A, D-H...A angle > 150 deg;
                acceptors: every chain O)
  hyd_core / hyd_glc   waters (O within 0.45 nm) per heavy atom of the imide core / glucose
Radial profiles: number density (atoms nm^-3) of core, linker, glucose, backbone heavy atoms of the largest aggregate, of pyrene heavy atoms
and of water O, vs distance from the largest aggregate's centre (0.1 nm shells), averaged over t >= 20 ns.
"""
import argparse, os, warnings
import numpy as np
import MDAnalysis as mda
from MDAnalysis.lib.distances import capped_distance, distance_array
from rdkit import Chem
warnings.filterwarnings('ignore')

ap = argparse.ArgumentParser()
ap.add_argument('--name', required=True); ap.add_argument('--gro', required=True); ap.add_argument('--xtc', default=None)
ap.add_argument('--npyr', type=int, required=True); ap.add_argument('--nchains', type=int, required=True)
ap.add_argument('--out', default='analysis/mechanism'); ap.add_argument('--params', default='params/monomer_c1')
ap.add_argument('--t0', type=float, default=20.0, help='start of the averaging window for profiles/contacts (ns)')
a = ap.parse_args()
xtc = a.xtc or f'site/work/{a.name}/whole500.xtc'
os.makedirs(a.out, exist_ok=True); NP, NC, P, C = a.npyr, a.nchains, 26, 556

# ---- groups from the monomer graph (same rules as render_many_movie.py)
mol = Chem.MolFromMolFile(os.path.join(a.params, 'monomer.mol'), removeHs=False)
raw = open(os.path.join(a.params, 'monomer_gaff2.mol2')).read().split('@<TRIPOS>ATOM')[1].split('@<TRIPOS>BOND')[0].strip().splitlines()
nm = [l.split()[1] for l in raw]
heavy_m = {i for i, at in enumerate(mol.GetAtoms()) if at.GetAtomicNum() > 1}
sp = mol.GetSubstructMatch(Chem.MolFromSmarts('[#7]-[CX4]-[CX4]-[CX4]-[CX4]-[OX2]')); linker = set(sp[1:])
five, pyran, imide = set(), None, None
for r in mol.GetRingInfo().AtomRings():
    if len(r) == 5:
        five |= set(r)
        if any(mol.GetAtomWithIdx(i).GetAtomicNum() == 7 for i in r): imide = [nm[i] for i in r]
    if len(r) == 6 and any(mol.GetAtomWithIdx(i).GetAtomicNum() == 8 for i in r): pyran = set(r)
glc = set(pyran); st = list(pyran)
while st:
    i = st.pop()
    for n in mol.GetAtomWithIdx(i).GetNeighbors():
        j = n.GetIdx()
        if j in heavy_m and j not in glc and j not in linker: glc.add(j); st.append(j)
vin = mol.GetSubstructMatches(Chem.MolFromSmarts('[CH2]=[CH1]-[CH1;R]')); backbone = {h[1] for h in vin} | {h[0] for h in vin}
core = (five | {n.GetIdx() for i in five for n in mol.GetAtomWithIdx(i).GetNeighbors() if n.GetAtomicNum() == 8 and n.GetIdx() not in five}) - backbone
GROUP = {nm[i]: ('backbone' if i in backbone else 'linker' if i in linker else 'glucose' if i in glc else 'core') for i in heavy_m}
GN = ['core', 'linker', 'glucose', 'backbone']

u = mda.Universe(a.gro, xtc)
names, resn, resid = u.atoms.names, u.atoms.resnames, u.atoms.resindices
nat = len(u.atoms); heavy = np.array([not n.startswith('H') for n in names])
mol_of = np.full(nat, -1)
for i in range(NP): mol_of[i * P:(i + 1) * P] = i
for k in range(NC): mol_of[NP * P + k * C: NP * P + (k + 1) * C] = NP + k
grp = np.full(nat, -1)
for i in range(NP * P, NP * P + NC * C):
    if heavy[i]: grp[i] = 3 if resn[i] in ('CPH', 'CPT') else GN.index(GROUP[names[i]])
ch_h = [np.where(heavy & (mol_of == NP + k))[0] for k in range(NC)]
py_h = [np.where(heavy & (mol_of == i))[0] for i in range(NP)]
all_ch_h = np.concatenate(ch_h) if NC else np.array([], int)
wat_o = np.where((resn == 'WAT') & (names == 'O'))[0]
# imide rings (one per monomer residue)
rings = []
for r in u.residues:
    if r.resname in ('MON', 'MNF', 'MNL'):
        idx = {n: i for n, i in zip(r.atoms.names, r.atoms.indices)}; rings.append([idx[n] for n in imide])
rings = np.array(rings, int).reshape(-1, 5); ring_chain = mol_of[rings[:, 0]] - NP if len(rings) else np.array([], int)
# glucose hydroxyl donors (O with a bonded H in the first frame) and acceptors (every chain O)
X0 = u.atoms.positions
ch_O = np.where((mol_of >= NP) & np.array([n.startswith('O') for n in names]))[0]
don_O, don_H = [], []
for o in ch_O:
    if grp[o] != 2: continue
    rr = np.where((resid == resid[o]) & np.array([n.startswith('H') for n in names]))[0]
    d = np.linalg.norm(X0[rr] - X0[o], axis=1)
    for hh in rr[d < 1.15]: don_O.append(o); don_H.append(hh)
don_O, don_H = np.array(don_O, int), np.array(don_H, int)
core_h = np.where(grp == 0)[0]; glc_h = np.where(grp == 2)[0]

def comps(adj):
    n = len(adj); lab = -np.ones(n, int); c = 0
    for s in range(n):
        if lab[s] >= 0: continue
        stq = [s]; lab[s] = c
        while stq:
            v = stq.pop()
            for w in np.where(adj[v] & (lab < 0))[0]: lab[w] = c; stq.append(w)
        c += 1
    return lab
def plane(Q):
    c = Q.mean(0); _, _, vt = np.linalg.svd(Q - c); return c, vt[2]
def mic(d, L): return d - L * np.round(d / L)
def stacked(cA, nA, cB, nB, L):
    d = mic(cB - cA, L); perp = abs(d @ nA); inpl = np.sqrt(max(d @ d - perp ** 2, 0))
    return abs(nA @ nB) > np.cos(np.radians(30)) and 3.0 <= perp <= 4.2 and inpl <= 4.0

edges = np.arange(0, 4.01, 0.1); shellV = 4 / 3 * np.pi * (edges[1:] ** 3 - edges[:-1] ** 3)
prof = {g: np.zeros(len(edges) - 1) for g in GN + ['pyrene', 'water']}; nprof = 0
contact_cnt = np.zeros(4); pairs_all = []; pp_all = []; agg_of_pyr = []
rows = []
for ts in u.trajectory:
    X = u.atoms.positions.copy(); L = ts.dimensions[:3].astype(float); box = ts.dimensions; t = ts.time / 1000
    r = dict(time_ns=round(t, 3))
    # molecule contact graph
    sol = np.where(heavy & (mol_of >= 0))[0]
    pr = capped_distance(X[sol], X[sol], 4.0, box=box, return_distances=False)
    mi, mj = mol_of[sol[pr[:, 0]]], mol_of[sol[pr[:, 1]]]; k_ = mi != mj
    A = np.zeros((NP + NC, NP + NC), bool); A[mi[k_], mj[k_]] = True; A |= A.T
    pairs_all.append(A[:NP, NP:].copy() if (NP and NC) else np.zeros((NP, NC), bool)); pp_all.append(A[:NP, :NP].copy())
    ag = np.full(NP, -1)
    if NC:
        lab = comps(A[NP:, NP:]); sizes = np.bincount(lab); big = int(np.argmax(sizes))
        # unwrap each aggregate: place chains in BFS order next to the running centroid
        com = np.array([X[h].mean(0) for h in ch_h]); agg_c, agg_R = {}, {}
        for c in range(lab.max() + 1):
            mem = list(np.where(lab == c)[0]); placed = [mem[0]]; cen = com[mem[0]].copy()
            for k in mem[1:]:
                sh = -L * np.round((com[k] - cen) / L); X[NP * P + k * C: NP * P + (k + 1) * C] += sh; com[k] += sh
                placed.append(k); cen = com[placed].mean(0)
            Q = np.concatenate([X[ch_h[k]] for k in mem]); cc = Q.mean(0); G = (Q - cc).T @ (Q - cc) / len(Q); ev = np.linalg.eigvalsh(G)
            agg_c[c] = cc; agg_R[c] = np.sqrt(5 / 3 * ev.sum())
            if c == big: r['agg_size'] = int(sizes[c]); r['agg_rg'] = np.sqrt(ev.sum()) / 10; r['agg_kappa2'] = 1.5 * (ev ** 2).sum() / ev.sum() ** 2 - 0.5
        r['n_aggregates'] = int(lab.max() + 1)
        # radial profiles around the largest aggregate (t >= t0)
        if t >= a.t0:
            cc = agg_c[big]; mem = np.where(lab == big)[0]
            for g in range(4):
                ids = np.concatenate([ch_h[k][grp[ch_h[k]] == g] for k in mem]); rd = np.linalg.norm(X[ids] - cc, axis=1) / 10
                prof[GN[g]] += np.histogram(rd, edges)[0]
            if NP: prof['pyrene'] += np.histogram(np.linalg.norm(mic(X[np.concatenate(py_h)] - cc, L), axis=1) / 10, edges)[0]
            near = capped_distance(cc[None], X[wat_o], 40.0, box=box, return_distances=True)
            prof['water'] += np.histogram(near[1] / 10, edges)[0]; nprof += 1
    # pyrene: depth, imide stacking, pyrene-pyrene stacking, contact chemistry
    if NP:
        pc = [plane(X[h]) for h in py_h]
        for i in range(NP):
            touch = np.where(A[i, NP:])[0]
            if NC and len(touch):
                c = np.bincount(lab[touch]).argmax(); ag[i] = c; r[f'depth_pyr{i+1}'] = np.linalg.norm(mic(pc[i][0] - agg_c[c], L)) / agg_R[c]
            else: r[f'depth_pyr{i+1}'] = np.nan
            ns = 0
            if len(rings):
                rc = X[rings].mean(1); cand = np.where(np.linalg.norm(mic(rc - pc[i][0], L), axis=1) < 6.0)[0]
                for q in cand:
                    cq, nq = plane(X[rings[q]])
                    ns += stacked(pc[i][0], pc[i][1], cq, nq, L)
            r[f'stack_imide_pyr{i+1}'] = ns
        agg_of_pyr.append(ag)
        r['pp_stacked'] = sum(stacked(pc[i][0], pc[i][1], pc[j][0], pc[j][1], L) for i in range(NP) for j in range(i))
        if NC and t >= a.t0:
            cp = capped_distance(X[np.concatenate(py_h)], X[all_ch_h], 4.0, box=box, return_distances=False)
            hit = np.unique(all_ch_h[cp[:, 1]]); contact_cnt += np.bincount(grp[hit], minlength=4)
    # hydrogen bonds between glucose OH and chain O
    if NC and len(don_O):
        pr = capped_distance(X[don_O], X[ch_O], 3.5, box=box, return_distances=False)
        dO, aO, dH = don_O[pr[:, 0]], ch_O[pr[:, 1]], don_H[pr[:, 0]]
        ok = dO != aO
        v1 = mic(X[dO] - X[dH], L); v2 = mic(X[aO] - X[dH], L)
        ang = np.degrees(np.arccos(np.clip((v1 * v2).sum(1) / np.linalg.norm(v1, axis=1) / np.linalg.norm(v2, axis=1), -1, 1)))
        ok &= ang > 150; inter = mol_of[dO] != mol_of[aO]
        r['hb_inter'] = int((ok & inter).sum()); r['hb_intra'] = int((ok & ~inter).sum())
        for nmx, ids in (('hyd_core', core_h), ('hyd_glc', glc_h)):
            pr2 = capped_distance(X[ids], X[wat_o], 4.5, box=box, return_distances=False); r[nmx] = len(pr2) / len(ids)
    rows.append(r)

import pandas as pd
df = pd.DataFrame(rows); df.to_csv(f'{a.out}/{a.name}_ts.csv', index=False, float_format='%.4f')
np.savez_compressed(f'{a.out}/{a.name}_pairs.npz', t=df.time_ns.values, pairs=np.array(pairs_all), pp=np.array(pp_all), agg_of_pyr=np.array(agg_of_pyr), patterns=np.array([f's{(k % 4) + 1}' for k in range(NC)]))
if NC and nprof:
    half = L.min() / 20; ok_r = edges[1:] <= half                     # beyond half the box edge the minimum-image shell is incomplete
    pd.DataFrame(dict(r_nm=(edges[:-1] + edges[1:]) / 2, **{g: np.where(ok_r, prof[g] / shellV / nprof, np.nan) for g in prof})).to_csv(f'{a.out}/{a.name}_profiles.csv', index=False, float_format='%.4f')
if NP and NC:
    tot = np.bincount(grp[all_ch_h], minlength=4)
    pd.DataFrame(dict(group=GN, chain_atoms=tot, atom_frac=tot / tot.sum(), contact_atoms=contact_cnt, contact_frac=contact_cnt / max(contact_cnt.sum(), 1),
                      enrichment=(contact_cnt / max(contact_cnt.sum(), 1)) / (tot / tot.sum()))).to_csv(f'{a.out}/{a.name}_contacts.csv', index=False, float_format='%.4f')
print(f'{a.name}: {len(df)} frames, columns {len(df.columns)}; donors {len(don_O)}; imide rings {len(rings)}')
