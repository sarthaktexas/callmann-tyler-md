import re
raw = open("monomer_gaff2.mol2").read()
sec = {}
for p in re.split(r"@<TRIPOS>", raw)[1:]:
    h, body = p.split("\n", 1)
    sec[h.strip()] = body.strip().split("\n")
atoms = [l.split() for l in sec["ATOM"]]
bonds = [l.split() for l in sec["BOND"]]
nm = {a[0]: a[1] for a in atoms}
info = {a[1]: dict(x=a[2], y=a[3], z=a[4], t=a[5], q=float(a[8])) for a in atoms}
bn = [(nm[b[1]], nm[b[2]], b[3]) for b in bonds]
q = {n: d["q"] for n, d in info.items()}

capatoms = ["C5", "H4", "H5", "C9", "H10", "H11"]
cap_each = sum(q[n] for n in capatoms) / 2
capC = (q["C5"] + q["C9"]) / 2
capH = (q["H4"] + q["H5"] + q["H10"] + q["H11"]) / 4

pairs = [("O","O7"),("C","C20"),("C1","C2"),("H","H1"),
         ("C3","C7"),("C4","C8"),("H2","H8"),("H3","H9")]
for a, b in pairs:
    q[a] = q[b] = (q[a] + q[b]) / 2
q["C4"] += cap_each
q["C8"] += cap_each
core = [a[1] for a in atoms if a[1] not in capatoms]
resid = -sum(q[n] for n in core)
q["C4"] += resid / 2
q["C8"] += resid / 2

def write(fname, res, alist, blist, head, tail):
    with open(fname, "w") as f:
        f.write(f"@<TRIPOS>MOLECULE\n{res}\n {len(alist)} {len(blist)} 1 0 0\nSMALL\nUSER_CHARGES\n\n@<TRIPOS>ATOM\n")
        for i, (n, x, y, z, t, c) in enumerate(alist, 1):
            f.write(f"{i:>6} {n:<5}{x:>10}{y:>10}{z:>10} {t:<4} 1 {res} {c:>13.8f}\n")
        f.write("@<TRIPOS>BOND\n")
        idx = {a[0]: i for i, a in enumerate(alist, 1)}
        for i, (a, b, o) in enumerate(blist, 1):
            f.write(f"{i:>6}{idx[a]:>6}{idx[b]:>6} {o}\n")
        f.write(f"@<TRIPOS>SUBSTRUCTURE\n     1 {res}         1 TEMP              0 ****  ****    0 ROOT\n")
        f.write(f"@<TRIPOS>HEADTAIL\n{head}\n{tail}\n")
    print(f"{fname}: {len(alist)} atoms, charge {sum(a[5] for a in alist):+.8f}")

def residue(res, qq):
    al = [(n, info[n]["x"], info[n]["y"], info[n]["z"], info[n]["t"], qq[n]) for n in core]
    bl = [b for b in bn if b[0] in core and b[1] in core]
    write(f"res_{res}.mol2", res, al, bl, "C4 1", "C8 1")

qF = dict(q); qF["C4"] -= cap_each
qL = dict(q); qL["C8"] -= cap_each
residue("MON", q); residue("MNF", qF); residue("MNL", qL)

def cap(res, c, h1, h2, head, tail):
    al = [("CC", info[c]["x"], info[c]["y"], info[c]["z"], info[c]["t"], capC),
          ("HC1", info[h1]["x"], info[h1]["y"], info[h1]["z"], info[h1]["t"], capH),
          ("HC2", info[h2]["x"], info[h2]["y"], info[h2]["z"], info[h2]["t"], capH)]
    write(f"res_{res}.mol2", res, al, [("CC","HC1","1"),("CC","HC2","1")], head, tail)

cap("CPH", "C5", "H4", "H5", "0 0", "CC 1")    # sits before the first monomer
cap("CPT", "C9", "H10", "H11", "CC 1", "0 0")  # sits after the last monomer
