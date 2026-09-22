#!/usr/bin/env python3
"""Generate ORCA vpot points and convert ORCA ESP output to Amber RESP ESP."""

import argparse
import math
from pathlib import Path

ANG_TO_BOHR = 1.8897259886

VDW = {
    "H": 1.20,
    "C": 1.70,
    "N": 1.55,
    "O": 1.52,
    "F": 1.47,
    "P": 1.80,
    "S": 1.80,
    "CL": 1.75,
    "BR": 1.85,
    "I": 1.98,
}


def read_xyz(path: Path):
    lines = path.read_text().splitlines()
    natom = int(lines[0].split()[0])
    atoms = []
    for line in lines[2 : 2 + natom]:
        fields = line.split()
        atoms.append((fields[0], tuple(float(x) for x in fields[1:4])))
    return atoms


def fibonacci_sphere(n: int):
    golden = math.pi * (3.0 - math.sqrt(5.0))
    for i in range(n):
        z = 1.0 - (2.0 * (i + 0.5) / n)
        r = math.sqrt(max(0.0, 1.0 - z * z))
        theta = golden * i
        yield (math.cos(theta) * r, math.sin(theta) * r, z)


def dist(a, b):
    return math.sqrt(sum((x - y) ** 2 for x, y in zip(a, b)))


def shell_points(atoms, density: float, scales):
    points = []
    for scale in scales:
        for idx, (elem, center) in enumerate(atoms):
            radius = VDW.get(elem.upper(), 1.70) * scale
            npts = max(50, int(4.0 * math.pi * radius * radius * density))
            for unit in fibonacci_sphere(npts):
                point = tuple(center[j] + radius * unit[j] for j in range(3))
                keep = True
                for j, (other_elem, other_center) in enumerate(atoms):
                    if j == idx:
                        continue
                    other_radius = VDW.get(other_elem.upper(), 1.70) * scale
                    if dist(point, other_center) < other_radius:
                        keep = False
                        break
                if keep:
                    points.append(point)
    return points


def write_orca_points(points, path: Path):
    with path.open("w") as handle:
        handle.write(f"{len(points)}\n")
        for x, y, z in points:
            handle.write(f"{x * ANG_TO_BOHR:16.8f}{y * ANG_TO_BOHR:16.8f}{z * ANG_TO_BOHR:16.8f}\n")


def read_vpot(path: Path):
    values = []
    for line in path.read_text().splitlines()[1:]:
        fields = line.split()
        if not fields:
            continue
        values.append(float(fields[-1]))
    return values


def write_resp_esp(atoms, points, potentials, path: Path):
    if len(points) != len(potentials):
        raise SystemExit(f"point/potential count mismatch: {len(points)} vs {len(potentials)}")
    with path.open("w") as handle:
        handle.write(f"{len(atoms):5d}{len(points):6d}\n")
        for _elem, xyz in atoms:
            handle.write("".join(f"{coord * ANG_TO_BOHR:16.7E}" for coord in xyz) + "\n")
        for xyz, pot in zip(points, potentials):
            handle.write(f"{pot:16.7E}" + "".join(f"{coord * ANG_TO_BOHR:16.7E}" for coord in xyz) + "\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--xyz", default="monomer_hf_opt.xyz")
    parser.add_argument("--points", default="monomer_resp.vpot.xyz")
    parser.add_argument("--vpot-out", default="monomer_resp.vpot.out")
    parser.add_argument("--resp-esp", default="monomer_resp.esp")
    parser.add_argument("--density", type=float, default=1.4, help="points per square Angstrom")
    parser.add_argument("--scales", default="1.4,1.6,1.8,2.0")
    parser.add_argument("--convert", action="store_true")
    args = parser.parse_args()

    atoms = read_xyz(Path(args.xyz))
    scales = [float(x) for x in args.scales.split(",") if x.strip()]
    points = shell_points(atoms, args.density, scales)

    if args.convert:
        potentials = read_vpot(Path(args.vpot_out))
        write_resp_esp(atoms, points, potentials, Path(args.resp_esp))
        print(f"wrote {args.resp_esp} with {len(atoms)} atoms and {len(points)} ESP points")
    else:
        write_orca_points(points, Path(args.points))
        print(f"wrote {args.points} with {len(points)} ESP points")


if __name__ == "__main__":
    main()
