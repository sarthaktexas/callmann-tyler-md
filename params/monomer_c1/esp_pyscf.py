#!/usr/bin/env python3
"""Replacement for the ORCA orca_vpot ESP step.

  esp   : HF/6-31G* single point (Cartesian d, as Gaussian/ORCA use) on the ORCA-optimised xyz,
          ESP on Merz-Kollman-style shells, written in Amber RESP .esp format.
  dipole: compare the dipole of a RESP charge list with the QM dipole (gate before using charges).

Why: the orca_vpot-derived pyrene .esp was ~1.96x the true HF/6-31G* potential (checked against
this very code), which inflated every RESP charge. Always run the dipole gate.
"""
import argparse, sys, numpy as np
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
from pyscf import gto, scf
from mk_resp_grid_from_orca import read_xyz, shell_points, ANG_TO_BOHR, write_resp_esp   # reuse grid code

A2B = ANG_TO_BOHR
def run_scf(xyz):
    atoms = read_xyz(Path(xyz))
    mol = gto.M(atom=[(e, x) for e, x in atoms], basis='6-31g*', unit='Angstrom', cart=True, verbose=0)
    mf = scf.RHF(mol).density_fit(); mf.conv_tol = 1e-9; mf.kernel()
    assert mf.converged, 'SCF did not converge'
    return atoms, mol, mf

def esp_at(mol, dm, pts):
    Z, C = mol.atom_charges(), mol.atom_coords()
    out = np.zeros(len(pts))
    for i in range(0, len(pts), 200):
        sl = slice(i, i + 200)
        out[sl] = -np.einsum('gij,ij->g', mol.intor('int1e_grids', grids=pts[sl]), dm)
    return out + sum(z / np.linalg.norm(pts - c, axis=1) for c, z in zip(C, Z))

def cmd_esp(a):
    atoms, mol, mf = run_scf(a.xyz)
    dm = mf.make_rdm1()
    S = mol.intor('int1e_ovlp')
    ne = float(np.trace(dm @ S)); nz = float(mol.atom_charges().sum())
    assert abs(ne - nz) < 1e-5, f'electron count {ne} != nuclear charge {nz}'
    pts_ang = shell_points(atoms, a.density, [float(x) for x in a.scales.split(',')])
    pts = np.array(pts_ang) * A2B
    V = esp_at(mol, dm, pts)
    far = esp_at(mol, dm, np.array([[100., 0, 0], [0, 100., 0]]))
    assert np.abs(far).max() < 1e-3, f'far-field potential {far} not ~0'
    write_resp_esp(atoms, pts_ang, V.tolist(), Path(a.out))
    print(f'E(HF/6-31G*) = {mf.e_tot:.8f}  points {len(pts)}  V range {V.min():.4f}..{V.max():.4f} a.u.')
    print(f'QM dipole (Debye) {np.round(mf.dip_moment(unit="Debye", verbose=0), 3)}  |mu| = {np.linalg.norm(mf.dip_moment(unit="Debye", verbose=0)):.3f}')
    print(f'wrote {a.out}')

def cmd_dipole(a):
    atoms, mol, mf = run_scf(a.xyz)
    mu_qm = np.array(mf.dip_moment(unit='Debye', verbose=0))
    q = np.array([float(x) for l in open(a.charges) for x in l.split() if not l.startswith('#')][:len(atoms)]) \
        if a.raw else np.loadtxt(a.charges, usecols=2, comments='#')
    X = np.array([x for _, x in atoms]); mu_q = (q[:, None] * X).sum(0) * 4.80320
    print('QM  dipole |mu| %.3f D  %s' % (np.linalg.norm(mu_qm), np.round(mu_qm, 3)))
    print('RESP dipole |mu| %.3f D  %s' % (np.linalg.norm(mu_q), np.round(mu_q, 3)))
    r = np.linalg.norm(mu_q) / max(np.linalg.norm(mu_qm), 1e-9)
    print(f'ratio {r:.2f}  ->', 'PASS' if 0.8 <= r <= 1.25 else 'FAIL (expected 0.8-1.25; ~2 means a 2x ESP)')
    sys.exit(0 if 0.8 <= r <= 1.25 else 1)

p = argparse.ArgumentParser(); s = p.add_subparsers(dest='cmd', required=True)
e = s.add_parser('esp'); e.add_argument('xyz'); e.add_argument('-o', '--out', default='monomer_resp.esp')
e.add_argument('--density', type=float, default=1.4); e.add_argument('--scales', default='1.4,1.6,1.8,2.0')
d = s.add_parser('dipole'); d.add_argument('xyz'); d.add_argument('charges'); d.add_argument('--raw', action='store_true')
a = p.parse_args(); {'esp': cmd_esp, 'dipole': cmd_dipole}[a.cmd](a)
