"""Pick one pull-trajectory frame per umbrella window and write the per-window mdp files.
Needs: pull_pullx.xvg, pull.xtc, pull.tpr (from 01_prep_pull.sh) and gmx on PATH. Standard library only.
Windows: 0.65-1.00 nm every 0.05 (starting distance is 0.655 nm), then 1.1-3.0 nm every 0.1 (k = 1000 kJ/mol/nm^2, sigma_x ~ 0.05 nm).
"""
import os, re, subprocess, sys
K = 1000.0
centres = [round(0.65 + 0.05 * i, 3) for i in range(8)] + [round(1.1 + 0.1 * i, 3) for i in range(20)]
os.makedirs('windows', exist_ok=True)
t, x = [], []
for l in open('pull_pullx.xvg'):
    if l[0] in '#@': continue
    f = l.split(); t.append(float(f[0])); x.append(float(f[1]))
print('pull distance range %.3f - %.3f nm over %d frames' % (min(x), max(x), len(x)))
gmx = os.environ.get('GMX', 'gmx')
tmpl = open('mdp/umbrella_template.mdp').read()
used = []
for i, d in enumerate(centres):
    j = min(range(len(x)), key=lambda k: abs(x[k] - d))
    off = x[j] - d
    if abs(off) > 0.05:
        sys.exit('window %d (d=%.2f): nearest pull frame is %.3f nm away - extend the pull' % (i, d, off))
    tag = '%02d' % i
    subprocess.run([gmx, 'trjconv', '-f', 'pull.xtc', '-s', 'pull.tpr', '-dump', str(t[j]), '-o', 'windows/conf_%s.gro' % tag],
                   input=b'0\n', check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    open('windows/umbrella_%s.mdp' % tag, 'w').write(
        tmpl.replace('@D@', '%.3f' % d).replace('@K@', '%g' % K).replace('@SEED@', str(1000 + i)))
    used.append('%s %.3f frame_t=%.1f ps start_d=%.3f' % (tag, d, t[j], x[j]))
open('windows.txt', 'w').write('\n'.join(used) + '\n')
print('\n'.join(used))
