import math, random, sys
inp, outp, ncis, seed = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
random.seed(seed)
lines = open(inp).read().splitlines()
isat = lambda l: l.startswith(("ATOM", "HETATM"))
res  = lambda l: int(l[22:26])
name = lambda l: l[12:16].strip()
xyz  = lambda l: [float(l[30:38]), float(l[38:46]), float(l[46:54])]
put  = lambda l, p: l[:30] + f"{p[0]:8.3f}{p[1]:8.3f}{p[2]:8.3f}" + l[54:]
cis = sorted(random.sample(range(2, 11), ncis))
for j in cis:
    a = next(xyz(l) for l in lines if isat(l) and res(l) == j and name(l) == "C8")
    b = next(xyz(l) for l in lines if isat(l) and res(l) == j + 1 and name(l) == "C4")
    ax = [y - x for x, y in zip(a, b)]
    n = math.sqrt(sum(c * c for c in ax)); ax = [c / n for c in ax]
    new = []
    for l in lines:
        if isat(l) and res(l) > j:
            p = xyz(l); v = [p[i] - a[i] for i in range(3)]
            d = sum(v[i] * ax[i] for i in range(3))
            l = put(l, [a[i] + 2 * d * ax[i] - v[i] for i in range(3)])
        new.append(l)
    lines = new
open(outp, "w").write("\n".join(lines) + "\n")
print(f"{outp}: cis junctions (residue j to j+1) = {cis}  ({ncis} of 9)")
