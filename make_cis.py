import math, sys
lines = open(sys.argv[1]).read().splitlines()
isat = lambda l: l.startswith(("ATOM", "HETATM"))
xyz  = lambda l: [float(l[30:38]), float(l[38:46]), float(l[46:54])]
name = lambda l: l[12:16].strip()
res  = lambda l: int(l[22:26])
atoms = [l for l in lines if isat(l)]
get = lambda n, r: next(xyz(l) for l in atoms if name(l) == n and res(l) == r)
a, b = get("C8", 1), get("C4", 2)
ax = [y - x for x, y in zip(a, b)]
n = math.sqrt(sum(c * c for c in ax)); ax = [c / n for c in ax]
def rot180(p):
    v = [p[i] - a[i] for i in range(3)]
    d = sum(v[i] * ax[i] for i in range(3))
    return [a[i] + 2 * d * ax[i] - v[i] for i in range(3)]
out = []
for l in lines:
    if isat(l) and res(l) == 2:
        x, y, z = rot180(xyz(l))
        l = l[:30] + f"{x:8.3f}{y:8.3f}{z:8.3f}" + l[54:]
    out.append(l)
open(sys.argv[2], "w").write("\n".join(out) + "\n")
