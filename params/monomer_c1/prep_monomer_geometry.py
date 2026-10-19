#!/usr/bin/env python3
"""Build or validate the starting geometry for a ROMP glycomonomer.

  build : stereo-defined SMILES -> lowest-MMFF-energy ETKDG conformer -> .xyz
  check : .xyz (any source: Rowan, ORCA, ...) -> formula, bonds, stereo, sugar identity

Needs RDKit >= 2022.09 (rdDetermineBonds, rdCIPLabeler).  Only the SMILES decides
what molecule you get; the geometry is only a starting point for the HF/6-31G* opt.
"""
import argparse, sys
from rdkit import Chem
from rdkit.Chem import AllChem, rdDetermineBonds, rdMolDescriptors, rdCIPLabeler

def mirror(smi):
    return smi.replace('@@', '!').replace('@', '@@').replace('!', '@')

# beta-D-glucopyranose (6R, 2R,3R,4S,5S); O-methyl placed where the alkyl chain goes
_BASE = {  # (position, anomer) -> SMILES of the D-series reference
    ('1-O-alkyl', 'beta'):  'OC[C@H]1O[C@@H](OC)[C@H](O)[C@@H](O)[C@@H]1O',
    ('1-O-alkyl', 'alpha'): 'OC[C@H]1O[C@H](OC)[C@H](O)[C@@H](O)[C@@H]1O',
    ('3-O-alkyl', 'beta'):  'OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](OC)[C@@H]1O',
    ('3-O-alkyl', 'alpha'): 'OC[C@H]1O[C@H](O)[C@H](O)[C@@H](OC)[C@@H]1O',
}

def sugar_identity(mol):
    """List every (series, position, anomer) reference that is a chirality-aware substructure."""
    heavy = Chem.RemoveHs(mol)
    hits = []
    for (pos, anom), smi in _BASE.items():
        for series, s in (('D', smi), ('L', mirror(smi))):
            if heavy.HasSubstructMatch(Chem.MolFromSmiles(s), useChirality=True):
                hits.append(f'{series}-gluco, {pos}, {anom}')
    return hits

# Ring-opened norbornene-dicarboximide cores (N-methyl stand-in). Both are achiral (meso).
# Only these two are reachable from a real norbornene: cis-fused imide, cis-1,3 vinyls.
# Verified against RDKit-embedded norbornenes (MMFF energy, C=O...alkene vs C=O...CH2 distance).
_CORES = {
    'endo': 'C=C[C@@H]1C[C@H](C=C)[C@@H]2C(=O)N(C)C(=O)[C@@H]21',
    'exo':  'C=C[C@@H]1C[C@H](C=C)[C@H]2C(=O)N(C)C(=O)[C@H]21',
}

def core_identity(mol):
    heavy = Chem.RemoveHs(mol)
    return [k for k, s in _CORES.items() if heavy.HasSubstructMatch(Chem.MolFromSmiles(s), useChirality=True)]

def report(mol, label):
    ok = True
    print(f'[{label}] formula {rdMolDescriptors.CalcMolFormula(mol)}, {mol.GetNumAtoms()} atoms')
    h = Chem.RemoveHs(mol)
    print(f'[{label}] SMILES (heavy atoms) {Chem.MolToSmiles(h)}')
    rdCIPLabeler.AssignCIPLabels(h)
    n_ster = sum(1 for a in h.GetAtoms() if a.HasProp('_CIPCode'))
    print(f'[{label}] CIP stereocentres: {n_ster}')
    ids = sugar_identity(mol)
    print(f'[{label}] sugar identity: {ids if ids else "NO D/L-gluco match - check the drawing"}')
    core = core_identity(mol)
    print(f'[{label}] norbornene-derived core: {core if core else "INVALID (not a cis-fused, cis-1,3-divinyl exo/endo core)"}')
    if not core:
        ok = False
    if not ids or any(i.startswith('L-') for i in ids):
        ok = False
        print(f'[{label}] WARNING: not D-glucose.')
    return ok

def cmd_build(a):
    m = Chem.AddHs(Chem.MolFromSmiles(a.smiles))
    cids = AllChem.EmbedMultipleConfs(m, a.nconf, randomSeed=a.seed)
    res = AllChem.MMFFOptimizeMoleculeConfs(m, maxIters=5000)
    best = min(range(len(res)), key=lambda i: res[i][1])
    conf = m.GetConformer(cids[best])
    Chem.AssignStereochemistryFrom3D(m, confId=cids[best])
    if Chem.MolToSmiles(Chem.RemoveHs(m)) != Chem.MolToSmiles(Chem.MolFromSmiles(a.smiles)):
        sys.exit('ERROR: stereochemistry changed during embedding')
    with open(a.out, 'w') as f:
        f.write(f'{m.GetNumAtoms()}\nsource: RDKit ETKDG+MMFF94 best of {a.nconf}; SMILES {a.smiles}\n')
        for at in m.GetAtoms():
            p = conf.GetAtomPosition(at.GetIdx())
            f.write(f'{at.GetSymbol():<2} {p.x:14.8f} {p.y:14.8f} {p.z:14.8f}\n')
    print(f'wrote {a.out}')
    report(m, 'built')

def cmd_mdl(a):
    """xyz (e.g. the ORCA HF/6-31G* result) -> MDL .mol with bond orders, for antechamber -fi mdl.
    Atom order is untouched, so it stays aligned with the RESP charge list."""
    L = open(a.xyz).read().splitlines(); n = int(L[0].split()[0])
    m = Chem.MolFromXYZBlock(f'{n}\n\n' + '\n'.join(L[2:2 + n]) + '\n')
    rdDetermineBonds.DetermineBonds(m, charge=0)
    Chem.AssignStereochemistryFrom3D(m)
    if a.smiles:
        if Chem.MolToSmiles(Chem.RemoveHs(m)) != Chem.MolToSmiles(Chem.MolFromSmiles(a.smiles)):
            sys.exit('ERROR: optimised geometry no longer matches the intended SMILES (stereo or connectivity)')
    Chem.MolToMolFile(m, a.out)
    print(f'wrote {a.out} ({n} atoms, atom order unchanged)')

def cmd_check(a):
    L = open(a.xyz).read().splitlines()
    n = int(L[0].split()[0])
    m = Chem.MolFromXYZBlock(f'{n}\n\n' + '\n'.join(L[2:2 + n]) + '\n')
    rdDetermineBonds.DetermineBonds(m, charge=a.charge)
    Chem.AssignStereochemistryFrom3D(m)
    ok = report(m, 'xyz')
    if a.smiles:
        same = Chem.MolToSmiles(Chem.RemoveHs(m)) == Chem.MolToSmiles(Chem.MolFromSmiles(a.smiles))
        print(f'[xyz] matches expected SMILES incl. stereo: {same}')
        ok = ok and same
    sys.exit(0 if ok else 1)

p = argparse.ArgumentParser()
s = p.add_subparsers(dest='cmd', required=True)
b = s.add_parser('build'); b.add_argument('smiles'); b.add_argument('-o', '--out', default='start.xyz')
b.add_argument('--nconf', type=int, default=300); b.add_argument('--seed', type=int, default=7)
m_ = s.add_parser('mdl'); m_.add_argument('xyz'); m_.add_argument('--smiles'); m_.add_argument('-o','--out',default='monomer.mol')
c = s.add_parser('check'); c.add_argument('xyz'); c.add_argument('--smiles'); c.add_argument('--charge', type=int, default=0)
a = p.parse_args()
{'build': cmd_build, 'check': cmd_check, 'mdl': cmd_mdl}[a.cmd](a)
