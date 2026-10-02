import re
raw = open("monomer_gaff2.mol2").read()
sec = {}
for p in re.split(r"@<TRIPOS>", raw)[1:]:
    h, body = p.split("\n", 1)
    sec[h.strip()] = body.strip().split("\n")
atoms = [l.split() for l in sec["ATOM"]]
bonds = [l.split() for l in sec["BOND"]]

q = {a[1]: float(a[8]) for a in atoms}
caps = ["C5", "H4", "H5", "C9", "H10", "H11"]
cap_each = sum(q[n] for n in caps) / 2

pairs = [("O","O7"), ("C","C20"), ("C1","C2"), ("H","H1"),
         ("C3","C7"), ("C4","C8"), ("H2","H8"), ("H3","H9")]
for a, b in pairs:
    q[a] = q[b] = (q[a] + q[b]) / 2
q["C4"] += cap_each
q["C8"] += cap_each

keep = [a for a in atoms if a[1] not in caps]
resid = -sum(q[a[1]] for a in keep)      # tiny rounding leftover
q["C4"] += resid / 2
q["C8"] += resid / 2

nid = {a[0]: i + 1 for i, a in enumerate(keep)}
kb = [b for b in bonds if b[1] in nid and b[2] in nid]
with open("monomer_resid.mol2", "w") as f:
    f.write(f"@<TRIPOS>MOLECULE\nMON\n {len(keep)} {len(kb)} 1 0 0\nSMALL\nUSER_CHARGES\n\n@<TRIPOS>ATOM\n")
    for a in keep:
        f.write(f"{nid[a[0]]:>6} {a[1]:<5}{a[2]:>10}{a[3]:>10}{a[4]:>10} {a[5]:<4} 1 MON {q[a[1]]:>10.6f}\n")
    f.write("@<TRIPOS>BOND\n")
    for i, b in enumerate(kb, 1):
        f.write(f"{i:>6}{nid[b[1]]:>6}{nid[b[2]]:>6} {b[3]}\n")
    f.write("@<TRIPOS>SUBSTRUCTURE\n     1 MON         1 TEMP              0 ****  ****    0 ROOT\n")
    f.write("@<TRIPOS>HEADTAIL\nC4 1\nC8 1\n")
print(len(keep), "atoms; charge sum =", sum(q[a[1]] for a in keep))
