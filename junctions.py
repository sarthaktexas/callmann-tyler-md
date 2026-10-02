import math, sys
at = {}
for l in open(sys.argv[1]):
    if l.startswith(("ATOM", "HETATM")):
        f = l.split()
        at[(f[2], int(f[4]))] = tuple(map(float, f[5:8]))
sub = lambda a, b: [x - y for x, y in zip(a, b)]
cr = lambda a, b: [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]
dot = lambda a, b: sum(x * y for x, y in zip(a, b))
def dih(p):
    b1, b2, b3 = sub(p[1], p[0]), sub(p[2], p[1]), sub(p[3], p[2])
    n1, n2 = cr(b1, b2), cr(b2, b3)
    m = cr(n1, [x / math.sqrt(dot(b2, b2)) for x in b2])
    return math.degrees(math.atan2(dot(m, n2), dot(n1, n2)))
for r in range(2, 11):
    p = [at[("C7", r)], at[("C8", r)], at[("C4", r + 1)], at[("C3", r + 1)]]
    d = sub(at[("C8", r)], at[("C4", r + 1)])
    print(f"junction {r}-{r+1}: dihedral {dih(p):7.1f}   C=C {math.sqrt(dot(d, d)):.3f} A")
