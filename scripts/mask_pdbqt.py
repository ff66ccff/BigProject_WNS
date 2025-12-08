#!/usr/bin/env python3
"""Mask receptor atoms in PDBQT file based on proximity to docked ligands."""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import List, Tuple


class AtomRecord:
    """Represents an atom in a PDBQT file."""
    
    def __init__(self, line: str) -> None:
        self.line = line.rstrip('\n')
        if len(line) >= 54:
            self.coord = (
                float(line[30:38]),
                float(line[38:46]),
                float(line[46:54]),
            )
        else:
            self.coord = (0.0, 0.0, 0.0)
        self.serial = int(line[6:11]) if line[6:11].strip() else 0
        self.is_atom = line.startswith(("ATOM", "HETATM"))
    
    def get_atom_type(self) -> str:
        """Get atom type from PDBQT line (columns 78-79)."""
        if len(self.line) >= 79:
            return self.line[77:79].strip()
        return ""
    
    def get_charge(self) -> float:
        """Get charge from PDBQT line (columns 71-76)."""
        if len(self.line) >= 77:
            try:
                return float(self.line[70:76].strip())
            except ValueError:
                return 0.0
        return 0.0
    
    def set_atom_type(self, atom_type: str) -> None:
        """Set atom type in PDBQT line (columns 78-79)."""
        if len(self.line) >= 79:
            self.line = self.line[:77] + atom_type.rjust(2) + self.line[79:]
        else:
            # Extend line if too short
            self.line = self.line.ljust(79)[:77] + atom_type.rjust(2)
    
    def set_charge(self, charge: float) -> None:
        """Set charge in PDBQT line (columns 71-76)."""
        charge_str = f"{charge:8.4f}"
        if len(self.line) >= 77:
            self.line = self.line[:70] + charge_str + self.line[78:]
        else:
            # Extend line if too short
            self.line = self.line.ljust(70) + charge_str


def calculate_distance(coord1: Tuple[float, float, float], 
                      coord2: Tuple[float, float, float]) -> float:
    """Calculate Euclidean distance between two coordinates."""
    return math.sqrt(
        (coord1[0] - coord2[0])**2 + 
        (coord1[1] - coord2[1])**2 + 
        (coord1[2] - coord2[2])**2
    )


def find_atoms_to_mask(receptor_atoms: List[AtomRecord], 
                      ligand_atoms: List[AtomRecord], 
                      cutoff: float = 3.5) -> List[int]:
    """Find receptor atoms within cutoff distance of any ligand atom."""
    atoms_to_mask = []
    
    for rec_idx, rec_atom in enumerate(receptor_atoms):
        if not rec_atom.is_atom:
            continue
            
        for lig_atom in ligand_atoms:
            if not lig_atom.is_atom:
                continue
                
            distance = calculate_distance(rec_atom.coord, lig_atom.coord)
            if distance < cutoff:
                atoms_to_mask.append(rec_idx)
                break  # No need to check other ligand atoms for this receptor atom
    
    return atoms_to_mask


def mask_receptor(receptor_file: Path, 
                  ligand_files: List[Path], 
                  output_file: Path, 
                  cutoff: float = 3.5) -> None:
    """Mask receptor atoms that are close to ligand atoms.
    
    Args:
        receptor_file: Path to receptor PDBQT file
        ligand_files: List of paths to ligand PDBQT files
        output_file: Path for output masked receptor PDBQT file
        cutoff: Distance cutoff in Angstroms for masking
    """
    # Load receptor atoms
    receptor_atoms = []
    for line in receptor_file.read_text(encoding="utf-8").splitlines():
        receptor_atoms.append(AtomRecord(line))
    
    # Load all ligand atoms
    all_ligand_atoms = []
    for ligand_file in ligand_files:
        for line in ligand_file.read_text(encoding="utf-8").splitlines():
            all_ligand_atoms.append(AtomRecord(line))
    
    # Find atoms to mask
    atoms_to_mask = find_atoms_to_mask(receptor_atoms, all_ligand_atoms, cutoff)
    
    # Write masked receptor
    with output_file.open('w', encoding="utf-8") as f:
        for idx, atom in enumerate(receptor_atoms):
            if idx in atoms_to_mask and atom.is_atom:
                # Mask this atom: set atom type to 'X' and charge to 0.000
                atom.set_atom_type('X')
                atom.set_charge(0.000)
            f.write(atom.line + '\n')
    
    print(f"Masked {len(atoms_to_mask)} receptor atoms in {output_file}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("receptor", type=Path, help="Receptor PDBQT file")
    parser.add_argument("ligands", nargs="+", type=Path, help="Ligand PDBQT files")
    parser.add_argument("-o", "--output", type=Path, required=True, 
                       help="Output masked receptor PDBQT file")
    parser.add_argument("-c", "--cutoff", type=float, default=3.5,
                       help="Distance cutoff in Angstroms (default: 3.5)")
    
    args = parser.parse_args()
    
    # Validate input files
    if not args.receptor.exists():
        raise FileNotFoundError(f"Receptor file not found: {args.receptor}")
    
    for ligand_file in args.ligands:
        if not ligand_file.exists():
            raise FileNotFoundError(f"Ligand file not found: {ligand_file}")
    
    mask_receptor(args.receptor, args.ligands, args.output, args.cutoff)


if __name__ == "__main__":
    main()