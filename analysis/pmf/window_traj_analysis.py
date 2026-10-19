"""Frame-level analysis of umbrella windows (solute-only xtc): pyrene coverage by polymer, chain Rg, contacts, per window.
coverage = 1 - SASA(pyrene in complex)/SASA(pyrene alone), solute only, Shrake-Rupley (300 pts, 0.14 nm probe, GROMACS-like radii);
chain made whole via bonds (fragment unwrap); pyrene shifted to the periodic image nearest the chain COG.
Usage: python window_traj_analysis.py <pattern_dir> <solute.tpr> <outdir> [stride=10] [wmax=17]   (stride 10 = every 100 ps)"""
import sys, os, numpy as np, MDAnalysis as mda
from concurrent.futures import ProcessPoolExecutor
P, TPR, OUT = sys.argv[1], sys.argv[2], sys.argv[3]; STRIDE = int(sys.argv[4]) if len(sys.argv) > 4 else 10; WMAX = int(sys.argv[5]) if len(sys.argv) > 5 else 17
RAD = {'C': 0.17, 'H': 0.11, 'N': 0.16, 'O': 0.15}; PROBE = 0.14
def golden(n):
    i = np.arange(n) + 0.5; phi = np.arccos(1 - 2 * i / n); th = np.pi * (1 + 5 ** 0.5) * i
    return np.c_[np.cos(th) * np.sin(phi), np.sin(th) * np.sin(phi), np.cos(phi)]
PTS = golden(300)
def sasa_pyr(X, r, pyr, occ):
    """X nm coords of all solute atoms; r radii+probe; pyr indices; occ indices of occluders (may include pyr)."""
    tot = 0.0
    for i in pyr:
        p = X[i] + PTS * r[i]; d = np.linalg.norm(p[:, None, :] - X[occ][None, :, :], axis=2)
        mask = occ != i; free = (d[:, mask] > r[occ][mask][None, :]).all(1); tot += 4 * np.pi * r[i] ** 2 * free.mean()
    return tot
def run(w):
    u = mda.Universe(TPR, f'{P}/windows/us_{w:02d}.xtc'); pyr = u.select_atoms('resname PYR'); chain = u.select_atoms('not resname PYR')
    el = np.array([n[0] for n in u.atoms.names]); r = np.array([RAD[e] for e in el]) + PROBE
    pi = pyr.indices; ci = chain.indices; heavy_c = chain.indices[np.array([n[0] != 'H' for n in chain.names])]; heavy_p = pyr.indices[np.array([n[0] != 'H' for n in pyr.names])]
    masses = chain.masses; rows = []
    for ts in u.trajectory[::STRIDE]:
        box = ts.dimensions[:3] / 10.0
        for f in u.atoms.fragments: f.unwrap(compound='fragments', reference='cog')
        pos = u.atoms.positions / 10.0
        c = pos[ci].mean(0); d = pos[pi].mean(0) - c; shift = -box * np.round(d / box); pos[pi] += shift
        near = np.where(np.linalg.norm(pos - pos[pi].mean(0), axis=2 if False else 1) < 1.0)[0]
        occ_all = np.union1d(near, pi); a_in = sasa_pyr(pos, r, pi, occ_all); a_alone = sasa_pyr(pos, r, pi, pi)
        com = (pos[ci] * masses[:, None]).sum(0) / masses.sum(); rg = np.sqrt((((pos[ci] - com) ** 2) * masses[:, None]).sum() / masses.sum())
        dm = np.linalg.norm(pos[heavy_c][:, None, :] - pos[heavy_p][None, :, :], axis=2).min(1)
        rows.append((ts.time / 1000.0, 1 - a_in / a_alone, rg, int((dm < 0.45).sum())))
    arr = np.array(rows); np.savetxt(f'{OUT}/window_{w:02d}.csv', arr, delimiter=',', header='t_ns,coverage,chain_Rg_nm,chain_heavy_within_0.45nm', comments='')
    return w, arr[:, 1].mean(), arr[:, 2].mean()
if __name__ == '__main__':
    os.makedirs(OUT, exist_ok=True)
    # one window per process (sandbox blocks multiprocessing): python window_traj_analysis.py ... <stride> <wmax> <single_window>
    ws = [int(sys.argv[6])] if len(sys.argv) > 6 else list(range(WMAX + 1))
    for w in ws:
        w, cv, rg = run(w); print('window %02d done: mean coverage %.2f, mean Rg %.2f nm' % (w, cv, rg), flush=True)
