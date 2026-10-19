#!/usr/bin/env python3
"""Build tleap repeat-unit residues for a ROMP norbornene-imide monomer from the molecular graph.

Replaces params/gaff2/make_residues.py, which hard-coded the first monomer's atom names.

  build_residues.py monomer.mol monomer_gaff2.mol2 [outdir]

monomer.mol        MDL file (prep_monomer_geometry.py mdl), vinyl-capped ring-opened monomer
monomer_gaff2.mol2 antechamber GAFF2 output (-c rc with the final RESP charges); SAME atom order as monomer.mol
Writes res_MON/MNF/MNL/CPH/CPT.mol2 in outdir and prints every check.

Convention (same as the first monomer): the two terminal =CH2 groups become the end caps CPH/CPT; the vinyl CH carbons become
HEAD/TAIL of the repeat unit. Steps: (1) average charges over topological symmetry classes (the two halves of the meso core);
(2) move half of the total cap charge onto each head/tail carbon; (3) tiny residual to make MON exactly neutral;
MNF/MNL = MON minus the cap charge on their head/tail carbon (the cap residue carries it).
"""
import sys, os, re, collections
from rdkit import Chem

molf, mol2f = sys.argv[1], sys.argv[2]
out = sys.argv[3] if len(sys.argv) > 3 else '.'
mol = Chem.MolFromMolFile(molf, removeHs=False)

raw = open(mol2f).read()
sec = {p.split('\n', 1)[0].strip(): p.split('\n', 1)[1].strip().split('\n') for p in re.split(r'@<TRIPOS>', raw)[1:]}
A = [l.split() for l in sec['ATOM']]
B = [l.split() for l in sec['BOND']]
assert len(A) == mol.GetNumAtoms(), 'atom count differs between .mol and .mol2'
for a, at in zip(A, mol.GetAtoms()):
    assert a[1][0].upper() == at.GetSymbol()[0].upper(), f'atom order mismatch at {a[1]}'
name = [a[1] for a in A]; typ = [a[5] for a in A]; xyz = [a[2:5] for a in A]; q = [float(a[8]) for a in A]
N = len(A)

# --- caps and heads from the graph
pat = Chem.MolFromSmarts('[CH2:1]=[CH1:2]-[CH1;R:3]')
hits = mol.GetSubstructMatches(pat)
assert len(hits) == 2, f'expected 2 terminal vinyl groups, found {len(hits)}'
capC = [h[0] for h in hits]; heads = [h[1] for h in hits]
capH = [n.GetIdx() for c in capC for n in mol.GetAtomWithIdx(c).GetNeighbors() if n.GetAtomicNum() == 1]
assert len(capH) == 4
capset = set(capC) | set(capH)
heads.sort(); head, tail = heads
cap_of = {head: [c for c in capC if mol.GetBondBetweenAtoms(c, head)][0], tail: [c for c in capC if mol.GetBondBetweenAtoms(c, tail)][0]}

# --- (1) symmetry averaging over topological classes
rank = list(Chem.CanonicalRankAtoms(mol, breakTies=False, includeChirality=False))
cls = collections.defaultdict(list)
for i, r in enumerate(rank): cls[r].append(i)
qs = list(q); nav = 0
for members in cls.values():
    if len(members) > 1:
        m = sum(q[i] for i in members) / len(members)
        for i in members: qs[i] = m
        nav += 1
sym_shift = max(abs(a - b) for a, b in zip(qs, q))

# --- (2) cap charge, (3) neutrality
cap_total = sum(qs[i] for i in capset)
cap_each = cap_total / 2
core = [i for i in range(N) if i not in capset]
qM = list(qs)
qM[head] += cap_each; qM[tail] += cap_each
resid = -sum(qM[i] for i in core)
qM[head] += resid / 2; qM[tail] += resid / 2
qF = list(qM); qF[head] -= cap_each
qL = list(qM); qL[tail] -= cap_each
capC_q = sum(qs[c] for c in capC) / 2
capH_q = sum(qs[h] for h in capH) / 4

def bonds_among(idx):
    s = set(idx); return [(int(b[1]) - 1, int(b[2]) - 1, b[3]) for b in B if int(b[1]) - 1 in s and int(b[2]) - 1 in s]

def write(fname, res, atoms, charges, bonds, headname, tailname):
    pos = {a: k + 1 for k, a in enumerate(atoms)}
    with open(os.path.join(out, fname), 'w') as f:
        f.write(f'@<TRIPOS>MOLECULE\n{res}\n {len(atoms)} {len(bonds)} 1 0 0\nSMALL\nUSER_CHARGES\n\n@<TRIPOS>ATOM\n')
        for a in atoms:
            f.write(f'{pos[a]:>6} {name[a]:<5}{xyz[a][0]:>10}{xyz[a][1]:>10}{xyz[a][2]:>10} {typ[a]:<4} 1 {res} {charges[a]:>13.8f}\n')
        f.write('@<TRIPOS>BOND\n')
        for k, (i, j, o) in enumerate(bonds, 1): f.write(f'{k:>6}{pos[i]:>6}{pos[j]:>6} {o}\n')
        f.write(f'@<TRIPOS>SUBSTRUCTURE\n     1 {res}         1 TEMP              0 ****  ****    0 ROOT\n@<TRIPOS>HEADTAIL\n{headname}\n{tailname}\n')
    print(f'{fname}: {len(atoms)} atoms, charge {sum(charges[a] for a in atoms):+.8f}')

hn, tn = f'{name[head]} 1', f'{name[tail]} 1'
for res, qq in (('MON', qM), ('MNF', qF), ('MNL', qL)):
    write(f'res_{res}.mol2', res, core, qq, bonds_among(core), hn, tn)

def cap(res, c, headname, tailname):
    hs = [n.GetIdx() for n in mol.GetAtomWithIdx(c).GetNeighbors() if n.GetAtomicNum() == 1]
    atoms = [c] + hs
    qq = {c: capC_q, hs[0]: capH_q, hs[1]: capH_q}
    names = {c: 'CC', hs[0]: 'HC1', hs[1]: 'HC2'}
    saved = list(name)
    for k, v in names.items(): name[k] = v
    write(f'res_{res}.mol2', res, atoms, qq, [(c, hs[0], '1'), (c, hs[1], '1')], headname, tailname)
    for k, v in enumerate(saved): name[k] = v
cap('CPH', cap_of[head], '0 0', 'CC 1'); cap('CPT', cap_of[tail], 'CC 1', '0 0')

print(f'head={name[head]} tail={name[tail]} cap carbons={[name[c] for c in capC]}')
print(f'symmetry classes averaged: {nav}, max charge change {sym_shift:.4f} e')
print(f'cap charge each {cap_each:+.5f} e, neutrality residual {resid:+.2e} e')
assert abs(sum(qM[i] for i in core)) < 1e-8
assert abs(sum(qF[i] for i in core) + cap_each) < 1e-8
# the cap + first residue must equal a neutral unit
assert abs((capC_q + 2 * capH_q) - cap_each) < 1e-8
print('checks passed: MON neutral, CPH + MNF neutral, CPT + MNL neutral')
