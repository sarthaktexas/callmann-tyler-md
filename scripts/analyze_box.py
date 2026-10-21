#!/usr/bin/env python3
"""Per-frame analysis of a pyrene + N chains box (or a chain-only box) with GROMACS tools.

  analyze_box.py --tpr md.tpr --xtc md.xtc --nchains 2 --out analysis/c1_2to1 [--dt 200] [--no-pyrene] [--chain-atoms 556]

Writes <out>.csv (per frame) and <out>_summary.txt. Columns: time_ns; pyrene_coverage_pct = 100*(1 - SASA_pyrene_in_complex / SASA_pyrene_alone) in the
SAME frame (0.14 nm probe, 960 points); mean chain Rg (nm) and per-chain Rg; minimum pyrene-chain atom distance (nm) per chain (gmx mindist, minimum image);
total chain SASA. Window means are given for 0-20 ns, 20 ns-end and the last 10 ns. State the starting state when quoting numbers.
Pitfalls handled: molecules made whole (-pbc whole), a hand-written index (no duplicate group names, ascending atom lists), SASA on a solute-only
trajectory with a solute-only .gro as structure (gmx sasa group selections and convert-tpr misbehaved on the 2:1 box otherwise), and no `mindist -pi`.
Solute atoms must be the first atoms in the system (pyrene = 1-26, then chains of --chain-atoms atoms each), which is how all boxes here are built.
Needs gmx on PATH (or GMX env var) and numpy.
"""
import argparse, os, subprocess, sys
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument('--tpr', required=True); ap.add_argument('--xtc', required=True); ap.add_argument('--nchains', type=int, required=True)
ap.add_argument('--out', required=True); ap.add_argument('--dt', type=float, default=200.0, help='frame stride in ps')
ap.add_argument('--no-pyrene', action='store_true'); ap.add_argument('--chain-atoms', type=int, default=556); ap.add_argument('--pyrene-atoms', type=int, default=26)
a = ap.parse_args()
GMX = os.environ.get('GMX', 'gmx'); out = os.path.abspath(a.out); work = out + '_work'; os.makedirs(work, exist_ok=True)
tpr, xtc = os.path.abspath(a.tpr), os.path.abspath(a.xtc)
npy = 0 if a.no_pyrene else a.pyrene_atoms; nsol = npy + a.nchains * a.chain_atoms
def gmx(args, stdin=None):
    r = subprocess.run(f'{GMX} {args}', shell=True, cwd=work, input=stdin, capture_output=True, text=True)
    return r
# index (full-system numbering == solute-only numbering because solute atoms come first)
groups = {}
if npy: groups['Pyrene'] = list(range(1, npy + 1))
groups['Solute'] = list(range(1, nsol + 1))
for k in range(a.nchains): groups[f'Chain{k + 1}'] = list(range(npy + k * a.chain_atoms + 1, npy + (k + 1) * a.chain_atoms + 1))
with open(f'{work}/idx.ndx', 'w') as f:
    for n, v in groups.items():
        f.write(f'[ {n} ]\n'); f.write('\n'.join(' '.join(map(str, v[i:i + 15])) for i in range(0, len(v), 15)) + '\n')
r = gmx(f'trjconv -s {tpr} -f {xtc} -pbc whole -dt {a.dt} -o whole.xtc', 'System\n'); assert os.path.exists(f'{work}/whole.xtc'), r.stderr[-400:]
r = gmx(f'trjconv -s {tpr} -f whole.xtc -n idx.ndx -o sol.xtc', 'Solute\n'); assert os.path.exists(f'{work}/sol.xtc'), r.stderr[-400:]
r = gmx(f'trjconv -s {tpr} -f whole.xtc -n idx.ndx -dump 0 -o sol.gro', 'Solute\n'); assert os.path.exists(f'{work}/sol.gro'), r.stderr[-400:]
load = lambda p: np.loadtxt(f'{work}/{p}', comments=('#', '@'), ndmin=2)
cols = {}
if npy:
    r = gmx('sasa -s sol.gro -f sol.xtc -n idx.ndx -surface \'all\' -output \'group "Pyrene"\' -probe 0.14 -ndots 960 -o sasa_complex.xvg'); assert os.path.exists(f'{work}/sasa_complex.xvg'), r.stderr[-400:]
    r = gmx('sasa -s sol.gro -f sol.xtc -n idx.ndx -surface \'group "Pyrene"\' -probe 0.14 -ndots 960 -o sasa_alone.xvg'); assert os.path.exists(f'{work}/sasa_alone.xvg'), r.stderr[-400:]
    c, al = load('sasa_complex.xvg'), load('sasa_alone.xvg'); t = c[:, 0]
    cols['pyrene_coverage_pct'] = 100 * (1 - c[:, 2] / al[:, 1]); cols['pyrene_sasa_alone_nm2'] = al[:, 1]
else:
    r = gmx('sasa -s sol.gro -f sol.xtc -n idx.ndx -surface \'all\' -probe 0.14 -ndots 960 -o sasa_chain.xvg'); c = load('sasa_chain.xvg'); t = c[:, 0]
if npy:
    r = gmx('sasa -s sol.gro -f sol.xtc -n idx.ndx -surface \'group "Solute"\' -output \'group "Solute"\' -probe 0.14 -ndots 960 -o sasa_solute.xvg')
rg = []
for k in range(a.nchains):
    r = gmx(f'gyrate -s {tpr} -f whole.xtc -n idx.ndx -o rg_{k}.xvg', f'Chain{k + 1}\n'); assert os.path.exists(f'{work}/rg_{k}.xvg'), r.stderr[-300:]
    rg.append(load(f'rg_{k}.xvg')[:, 1]); cols[f'Rg_chain{k + 1}_nm'] = rg[-1]
cols['Rg_mean_nm'] = np.mean(rg, axis=0)
if npy:
    for k in range(a.nchains):
        r = gmx(f'mindist -s {tpr} -f whole.xtc -n idx.ndx -od md_{k}.xvg', f'Pyrene\nChain{k + 1}\n'); assert os.path.exists(f'{work}/md_{k}.xvg'), r.stderr[-300:]
        cols[f'mindist_pyrene_chain{k + 1}_nm'] = load(f'md_{k}.xvg')[:, 1]
n = min(len(t), *[len(v) for v in cols.values()]); t = t[:n]
with open(out + '.csv', 'w') as f:
    f.write('time_ns,' + ','.join(cols) + '\n')
    for i in range(n): f.write(f'{t[i] / 1000:.3f},' + ','.join(f'{cols[k][i]:.4f}' for k in cols) + '\n')
tend = t[-1]; wins = [('0-20 ns', 0, 20000), ('20 ns-end', 20000, tend + 1), ('last 10 ns', tend - 10000, tend + 1), ('all', 0, tend + 1)]
lines = [f'{os.path.basename(out)}: {n} frames, stride {a.dt:.0f} ps, {t[0] / 1000:.1f}-{tend / 1000:.1f} ns; chains {a.nchains}; pyrene {"no" if a.no_pyrene else "yes"}']
for lab, lo, hi in wins:
    s = (t >= lo) & (t < hi)
    row = f'  {lab:11s}'
    if npy: row += f' coverage {cols["pyrene_coverage_pct"][s].mean():5.1f}% (min {cols["pyrene_coverage_pct"][s].min():5.1f}, max {cols["pyrene_coverage_pct"][s].max():5.1f})'
    row += f'  Rg {cols["Rg_mean_nm"][s].mean():.2f} nm'
    if npy: row += '  min dist ' + ', '.join(f'{cols[f"mindist_pyrene_chain{k + 1}_nm"][s].mean():.2f}' for k in range(a.nchains)) + ' nm'
    lines.append(row)
if npy: lines.append(f'  first frame coverage {cols["pyrene_coverage_pct"][0]:.1f}%  last frame {cols["pyrene_coverage_pct"][-1]:.1f}%')
open(out + '_summary.txt', 'w').write('\n'.join(lines) + '\n'); print('\n'.join(lines))
