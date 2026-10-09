"""Build the starting box for the pyrene-chain umbrella-sampling pilot (chain pattern 1).

Input : the final frame (100 ns) of the 1:1 seed-1 run (bound state), plus its tpr for bonds/elements.
Output: solute.gro  (pyrene + chain, made whole, centred in a cubic box)
        index.ndx   (System, PYR, Chain, Chain_core, Solute, Water_and_ions are written after solvation by make_index.py)
Usage : python build_start.py <md.tpr> <md.xtc> [box_nm=7.0]
"""
import sys, numpy as np, MDAnalysis as mda
from MDAnalysis.transformations import unwrap

tpr, xtc = sys.argv[1], sys.argv[2]
box = float(sys.argv[3]) if len(sys.argv) > 3 else 7.0
u = mda.Universe(tpr, xtc)
u.trajectory[-1]
print('frame time (ps):', u.trajectory.ts.time)
sol = u.select_atoms('not resname WAT')
assert sol.n_atoms == 582, sol.n_atoms
# make each molecule whole via the bond graph, then put pyrene in the periodic image nearest the chain
for frag in sol.fragments:
    frag.unwrap(compound='fragments', reference='cog')
pyr = sol.select_atoms('resname PYR'); chain = sol.select_atoms('not resname PYR')
box0 = u.dimensions[:3]
shift = np.round((chain.center_of_geometry() - pyr.center_of_geometry()) / box0) * box0
pyr.positions += shift
com = sol.center_of_geometry()
sol.positions += (np.array([box, box, box]) * 10 / 2 - com)
sol.write('solute_unboxed.gro')
print('pyrene-chain COG distance (nm):', np.linalg.norm(pyr.center_of_geometry() - chain.center_of_geometry()) / 10)
