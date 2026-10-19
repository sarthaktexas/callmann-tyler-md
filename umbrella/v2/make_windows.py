"""Pick one pull frame per window and write per-window mdp files. Run from the pattern directory.
Usage: python make_windows.py <in|out> <pull deffnm: pull_in|pull_out>
Windows: 0.55-1.00 nm every 0.05 and 1.1-1.8 every 0.1 (NEAR: 80 ns, first 30 ns = pre-equilibration, discarded),
         1.9-3.0 every 0.1 (FAR: 12 ns, first 2 ns discarded). k = 1000 kJ/mol/nm^2.
Writes windows_<dir>/conf_XX.gro, umbrella_XX.mdp and windows_<dir>.txt (idx centre frame_t start_d nsteps equil_ns)."""
import os, subprocess, sys
DIR, DEFFNM = sys.argv[1], sys.argv[2]; K = 1000.0
HERE = os.path.dirname(os.path.abspath(__file__)); PATTERN = int(os.environ['PATTERN'])
near = [round(0.55 + 0.05 * i, 3) for i in range(10)] + [round(1.1 + 0.1 * i, 3) for i in range(8)]
far = [round(1.9 + 0.1 * i, 3) for i in range(12)]
centres = near + far
NSTEPS_NEAR, EQ_NEAR, NSTEPS_FAR, EQ_FAR = 40_000_000, 30.0, 6_000_000, 2.0
wd = 'windows_' + DIR; os.makedirs(wd, exist_ok=True)
t, x = [], []
for l in open(DEFFNM + '_pullx.xvg'):
    if l[0] in '#@': continue
    f = l.split(); t.append(float(f[0])); x.append(float(f[1]))
print('%s: pull distance range %.3f - %.3f nm over %d frames' % (DEFFNM, min(x), max(x), len(x)))
gmx = os.environ.get('GMX', 'gmx'); tmpl = open(os.path.join(HERE, 'mdp/umbrella_template.mdp')).read()
rows = []
for i, d in enumerate(centres):
    j = min(range(len(x)), key=lambda k: abs(x[k] - d)); off = x[j] - d
    if abs(off) > 0.03: sys.exit('window %d (d=%.2f): nearest pull frame is %.3f nm away - lengthen/extend the pull' % (i, d, off))
    ns, eq = (NSTEPS_NEAR, EQ_NEAR) if d < 1.85 else (NSTEPS_FAR, EQ_FAR)
    tag = '%02d' % i; seed = 1000 * PATTERN + (0 if DIR == 'in' else 500) + i
    subprocess.run([gmx, 'trjconv', '-f', DEFFNM + '.xtc', '-s', DEFFNM + '.tpr', '-dump', str(t[j]), '-o', '%s/conf_%s.gro' % (wd, tag)],
                   input=b'0\n', check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    open('%s/umbrella_%s.mdp' % (wd, tag), 'w').write(tmpl.replace('@D@', '%.3f' % d).replace('@K@', '%g' % K).replace('@SEED@', str(seed)).replace('@NSTEPS@', str(ns)))
    rows.append('%s %.3f frame_t=%.1f ps start_d=%.3f nsteps=%d equil_ns=%g' % (tag, d, t[j], x[j], ns, eq))
open('windows_%s.txt' % DIR, 'w').write('\n'.join(rows) + '\n'); print('\n'.join(rows))
