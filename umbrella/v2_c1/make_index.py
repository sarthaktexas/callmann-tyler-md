"""Write index.ndx for a solvated pyrene+chain box: System, PYRENE (1-26), CHAIN (27-582), CORE (120), SOLUTE (1-582), Water_and_ions.
Usage: python make_index.py <sys.gro> <out index.ndx>"""
import sys, os
gro, out = sys.argv[1], sys.argv[2]
HERE = os.path.dirname(os.path.abspath(__file__))
n = int(open(gro).read().split('\n')[1])
core = [int(x) for x in open(os.path.join(HERE, 'core_heavy_atoms.txt')).read().split()]
def grp(name, ids):
    s = '[ %s ]\n' % name
    for i in range(0, len(ids), 15): s += ' '.join('%d' % v for v in ids[i:i + 15]) + '\n'
    return s
open(out, 'w').write(grp('System', list(range(1, n + 1))) + grp('PYRENE', list(range(1, 27))) + grp('CHAIN', list(range(27, 583)))
                     + grp('CORE', sorted(core)) + grp('SOLUTE', list(range(1, 583))) + grp('Water_and_ions', list(range(583, n + 1))))
print('index written: %d atoms, %d core' % (n, len(core)))
