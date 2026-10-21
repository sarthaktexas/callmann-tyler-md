"""Measure the CORE-PYRENE COM distance in a whole-molecule .gro and print the nsteps a slow pull needs to reach a target.
Usage: python plan_pull.py <whole.gro> <index.ndx> <target_nm> <rate_nm_per_ps> [dt_ps=0.002] [margin_nm=0.0]
Prints:  D0=<nm> NSTEPS=<int>     (use `eval` or parse). The gro must have whole molecules (gmx trjconv -pbc whole)."""
import sys, numpy as np, math
gro, ndx, target, rate = sys.argv[1], sys.argv[2], float(sys.argv[3]), float(sys.argv[4])
dt = float(sys.argv[5]) if len(sys.argv) > 5 else 0.002
L = open(gro).read().split('\n'); n = int(L[1]); box = np.array([float(v) for v in L[2 + n].split()[:3]])
names = [l[10:15].strip() for l in L[2:2 + n]]; xyz = np.array([[float(l[20:28]), float(l[28:36]), float(l[36:44])] for l in L[2:2 + n]])
MASS = {'H': 1.008, 'C': 12.011, 'N': 14.007, 'O': 15.999}
groups = {}; cur = None
for l in open(ndx):
    l = l.strip()
    if l.startswith('['): cur = l.strip('[] ').strip(); groups[cur] = []
    elif l and cur: groups[cur] += [int(v) - 1 for v in l.split()]
def com(ids):
    m = np.array([MASS.get(names[i][0], 12.0) for i in ids]); return (xyz[ids] * m[:, None]).sum(0) / m.sum()
d = com(groups['PYRENE']) - com(groups['CORE']); d -= box * np.round(d / box)
D0 = float(np.linalg.norm(d)); span = abs(target - D0)
nsteps = math.ceil(span / abs(rate) / dt)
print('D0=%.4f NSTEPS=%d' % (D0, nsteps))
