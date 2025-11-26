#!/usr/bin/env python3
"""Combine protein and ligand coordinates into a filtered complex."""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Iterable, List, Tuple


Coordinate = Tuple[float, float, float]


class AtomRecord:
    __slots__ = ("line", "coord", "serial")

    def __init__(self, line: str) -> None:
        self.line = line
        self.coord = (
            float(line[30:38]),
            float(line[38:46]),
            float(line[46:54]),
        )
        self.serial = int(line[6:11]) if line[6:11].strip() else 0


def load_atoms(path: Path) -> List[AtomRecord]:
    atoms: List[AtomRecord] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(("ATOM", "HETATM")):
            atoms.append(AtomRecord(line.ljust(80)))
    return atoms


def min_distance(atoms_a: Iterable[AtomRecord], atoms_b: Iterable[AtomRecord]) -> float:
    min_dist = float("inf")
    for atom_a in atoms_a:
        ax, ay, az = atom_a.coord
        for atom_b in atoms_b:
            bx, by, bz = atom_b.coord
            dist = math.sqrt((ax - bx) ** 2 + (ay - by) ** 2 + (az - bz) ** 2)
            if dist < min_dist:
                min_dist = dist
    return min_dist if min_dist != float("inf") else 0.0


def filter_ligands(ligand_atoms: List[List[AtomRecord]], cutoff: float) -> List[List[AtomRecord]]:
    filtered: List[List[AtomRecord]] = []
    for atoms in ligand_atoms:
        keep = True
        for accepted in filtered:
            if min_distance(atoms, accepted) < cutoff:
                keep = False
                break
        if keep:
            filtered.append(atoms)
    return filtered


def write_complex(
    protein_atoms: List[AtomRecord],
    ligand_groups: List[List[AtomRecord]],
    output_path: Path,
) -> None:
    serial = 1
    lines: List[str] = []
    for atom in protein_atoms:
        line = f"{atom.line[:6]}{serial:5d}{atom.line[11:]}"
        lines.append(line)
        serial += 1
    for group in ligand_groups:
        for atom in group:
            line = f"{atom.line[:6]}{serial:5d}{atom.line[11:]}"
            lines.append(line)
            serial += 1
    lines.append("END")
    output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("protein", type=Path, help="Path to protein PDB/PDBQT")
    parser.add_argument("ligands", nargs="+", type=Path, help="Ligand pose files (PDBQT)")
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path("../complex/complex_filtered.pdb"),
        help="Output path for merged complex",
    )
    parser.add_argument(
        "--cutoff",
        type=float,
        default=2.0,
        help="Distance cutoff in Å for removing clashing ligand poses",
    )
    args = parser.parse_args()

    protein_atoms = load_atoms(args.protein)
    ligand_atoms = [load_atoms(path) for path in args.ligands]
    ligand_groups = filter_ligands(ligand_atoms, args.cutoff)

    write_complex(protein_atoms, ligand_groups, args.output)


if __name__ == "__main__":
    main()
