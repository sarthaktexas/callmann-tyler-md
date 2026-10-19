#!/usr/bin/env python3
"""Write a copy of a GROMACS .top with the PYR [ atoms ] charges replaced by a new RESP set.

usage: swap_pyrene_charges.py IN.top OUT.top [NEW_CHARGES.dat] [OLD_CHARGES.dat]
Defaults: new = params/pyrene/pyrene_resp_charges_corrected.dat, old = params/pyrene/resp_charges_stage2.dat
Safety checks: atom count 26, atom order identical to the xyz-ordered charge files (old charges in the .top must
equal the old RESP list to 1e-6), the new charges are shifted by sum/26 so PYR is exactly neutral. Nothing else is touched.
"""
import re, sys, os
here = os.path.dirname(os.path.abspath(__file__))
P = os.path.join(here, '..', 'params', 'pyrene')
src, dst = sys.argv[1], sys.argv[2]
new_f = sys.argv[3] if len(sys.argv) > 3 else os.path.join(P, 'pyrene_resp_charges_corrected.dat')
old_f = sys.argv[4] if len(sys.argv) > 4 else os.path.join(P, 'resp_charges_stage2.dat')
rd = lambda f: [float(l.split()[2]) for l in open(f) if l.strip() and not l.startswith('#')]
new, old = rd(new_f), rd(old_f)
assert len(new) == len(old) == 26
shift = sum(new) / 26.0
new = [q - shift for q in new]
lines = open(src).read().split('\n')
i = next(k for k, l in enumerate(lines) if l.strip() == '[ moleculetype ]' and lines[k + 2].split()[0] == 'PYR')
a0 = next(k for k in range(i, len(lines)) if lines[k].strip() == '[ atoms ]') + 1
idx, out, qtot = 0, list(lines), 0.0
for k in range(a0, len(lines)):
    l = lines[k]
    if l.startswith('['):
        break
    if not l.strip() or l.lstrip().startswith(';'):
        continue
    f = l.split(';')[0].split()
    q_old = float(f[6])
    assert abs(q_old - old[idx]) < 1e-6, f'atom {idx+1}: top has {q_old}, old RESP list has {old[idx]} (different atom order?)'
    qtot += new[idx]
    body = re.sub(re.escape(f[6]), f'{new[idx]: .8f}'.strip(), l.split(';')[0], count=1)
    out[k] = body + f' ; qtot {qtot:.6f}'
    idx += 1
assert idx == 26, idx
open(dst, 'w').write('\n'.join(out))
print(f'{dst}: PYR charges replaced ({idx} atoms); H charges now', [round(new[j], 3) for j in range(16, 26)], 'sum', round(sum(new), 9))
