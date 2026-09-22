#!/usr/bin/env python3
"""Create two-stage Amber RESP input files from an XYZ geometry."""

import argparse
import math
from pathlib import Path

ATNO = {
    "H": 1,
    "C": 6,
    "N": 7,
    "O": 8,
    "F": 9,
    "P": 15,
    "S": 16,
    "CL": 17,
    "BR": 35,
    "I": 53,
}

COV = {
    "H": 0.31,
    "C": 0.76,
    "N": 0.71,
    "O": 0.66,
    "F": 0.57,
    "P": 1.07,
    "S": 1.05,
    "CL": 1.02,
    "BR": 1.20,
    "I": 1.39,
}


def read_xyz(path):
    lines = Path(path).read_text().splitlines()
    natom = int(lines[0].split()[0])
    atoms = []
    for line in lines[2 : 2 + natom]:
        fields = line.split()
        atoms.append((fields[0], tuple(float(x) for x in fields[1:4])))
    return atoms


def distance(a, b):
    return math.sqrt(sum((x - y) ** 2 for x, y in zip(a, b)))


def bonded(atoms, i, j):
    ei, xi = atoms[i]
    ej, xj = atoms[j]
    cutoff = COV.get(ei.upper(), 0.76) + COV.get(ej.upper(), 0.76) + 0.45
    return distance(xi, xj) <= cutoff


def h_equivalences(atoms):
    """Return RESP ivary list with hydrogens on the same heavy atom constrained."""
    ivary = [0] * len(atoms)
    attached = {}
    for i, (elem, _xyz) in enumerate(atoms):
        if elem.upper() != "H":
            continue
        heavy = None
        best = 999.0
        for j, (other, _other_xyz) in enumerate(atoms):
            if other.upper() == "H":
                continue
            if bonded(atoms, i, j):
                d = distance(atoms[i][1], atoms[j][1])
                if d < best:
                    best = d
                    heavy = j
        if heavy is not None:
            attached.setdefault(heavy, []).append(i)

    for hydrogens in attached.values():
        if len(hydrogens) < 2:
            continue
        reference = hydrogens[0] + 1
        for hidx in hydrogens[1:]:
            ivary[hidx] = reference
    return ivary


def write_resp_input(path, title, qwt, atoms, ivary, iqopt):
    with Path(path).open("w") as handle:
        handle.write(f"{title}\n")
        handle.write(" &cntrl\n")
        handle.write("  nmol = 1,\n")
        handle.write("  ihfree = 1,\n")
        handle.write(f"  iqopt = {iqopt},\n")
        handle.write(f"  qwt = {qwt:.6f}\n")
        handle.write(" &end\n")
        handle.write("  1.0\n")
        handle.write("Monomer\n")
        handle.write(f"{0:5d}{len(atoms):5d}\n")
        for (elem, _xyz), iv in zip(atoms, ivary):
            handle.write(f"{ATNO[elem.upper()]:5d}{iv:5d}\n")
        handle.write("\n")


def write_qin(path, atoms, charges=None):
    if charges is None:
        charges = [0.0] * len(atoms)
    with Path(path).open("w") as handle:
        for i in range(0, len(charges), 8):
            handle.write("".join(f"{q:10.6f}" for q in charges[i : i + 8]) + "\n")


def read_charges(path, natom):
    path = Path(path)
    if not path.exists() and path.suffix == ".qout":
        out_path = path.with_suffix(".out")
        if out_path.exists():
            return read_resp_out_charges(out_path, natom)
    vals = []
    for field in path.read_text().split():
        try:
            vals.append(float(field))
        except ValueError:
            pass
    if len(vals) < natom:
        out_path = path.with_suffix(".out")
        if path.suffix == ".qout" and out_path.exists():
            return read_resp_out_charges(out_path, natom)
        raise SystemExit(f"Could not read {natom} charges from {path}")
    return vals[:natom]


def read_resp_out_charges(path, natom):
    charges = []
    in_table = False
    for line in Path(path).read_text(errors="ignore").splitlines():
        if "q(opt)" in line and "ivary" in line:
            in_table = True
            charges = []
            continue
        if not in_table:
            continue
        fields = line.split()
        if len(fields) < 4:
            if charges:
                break
            continue
        try:
            int(fields[0])
            charge = float(fields[2])
        except ValueError:
            continue
        charges.append(charge)
        if len(charges) == natom:
            return charges
    if len(charges) < natom:
        raise SystemExit(f"Could not read {natom} q(opt) charges from {path}")
    return charges[:natom]


def write_charge_table(path, atoms, charges):
    with Path(path).open("w") as handle:
        handle.write("# idx element charge\n")
        for i, ((elem, _xyz), charge) in enumerate(zip(atoms, charges), 1):
            handle.write(f"{i:5d} {elem:2s} {charge: .8f}\n")
        handle.write(f"# charge_sum {sum(charges): .10f}\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--xyz", default="monomer_hf_opt.xyz")
    parser.add_argument("--stage1-qout")
    parser.add_argument("--stage2-qout")
    args = parser.parse_args()

    atoms = read_xyz(args.xyz)
    if args.stage2_qout:
        charges = read_charges(args.stage2_qout, len(atoms))
        write_charge_table("resp_charges_stage2.dat", atoms, charges)
        print(f"wrote resp_charges_stage2.dat; charge sum = {sum(charges):.10f}")
        return

    stage1_ivary = [0] * len(atoms)
    stage2_ivary = h_equivalences(atoms)
    write_resp_input("resp_stage1.in", "Stage 1 RESP fit", 0.0005, atoms, stage1_ivary, 1)
    write_qin("resp_stage1.qin", atoms)

    if args.stage1_qout:
        stage1 = read_charges(args.stage1_qout, len(atoms))
        write_resp_input("resp_stage2.in", "Stage 2 RESP fit", 0.0010, atoms, stage2_ivary, 2)
        write_qin("resp_stage2.qin", atoms, stage1)
    else:
        write_resp_input("resp_stage2.in", "Stage 2 RESP fit", 0.0010, atoms, stage2_ivary, 2)
        write_qin("resp_stage2.qin", atoms)

    constrained = sum(1 for x in stage2_ivary if x != 0)
    print(f"wrote RESP inputs for {len(atoms)} atoms; stage-2 constrained atoms: {constrained}")


if __name__ == "__main__":
    main()
