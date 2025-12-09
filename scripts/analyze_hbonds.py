#!/usr/bin/env python3
"""Hydrogen bond analysis for WnS project: Count H-bonds and filter stable interactions."""

from __future__ import annotations

import argparse
import csv
import math
import os
import sys
from pathlib import Path
from typing import Dict, List, Tuple, Set, Optional

try:
    import numpy as np
except ImportError:
    print("Warning: numpy not available, using slower pure Python calculations")
    np = None

try:
    import yaml  # type: ignore
except ImportError as exc:  # pragma: no cover
    raise SystemExit("PyYAML is required. Install it via 'pip install pyyaml'.") from exc


class Atom:
    """Represents an atom in a PDB file."""
    
    def __init__(self, line: str) -> None:
        self.line = line.rstrip('\n')
        self.atom_id = int(line[6:11])
        self.atom_name = line[12:16].strip()
        self.residue_name = line[17:20].strip()
        self.residue_id = int(line[22:26])
        self.coord = (
            float(line[30:38]),
            float(line[38:46]),
            float(line[46:54]),
        )
        self.element = line[76:78].strip() if len(line) > 76 else self.atom_name[0]


class HydrogenBond:
    """Represents a hydrogen bond between donor and acceptor."""
    
    def __init__(self, donor: Atom, acceptor: Atom, distance: float, angle: float) -> None:
        self.donor = donor
        self.acceptor = acceptor
        self.distance = distance
        self.angle = angle


def calculate_distance(coord1: Tuple[float, float, float], 
                      coord2: Tuple[float, float, float]) -> float:
    """Calculate Euclidean distance between two coordinates."""
    if np:
        return np.linalg.norm(np.array(coord1) - np.array(coord2))
    else:
        return math.sqrt(
            (coord1[0] - coord2[0])**2 + 
            (coord1[1] - coord2[1])**2 + 
            (coord1[2] - coord2[2])**2
        )


def calculate_angle(donor_h: Tuple[float, float, float], 
                   donor: Tuple[float, float, float],
                   acceptor: Tuple[float, float, float]) -> float:
    """Calculate angle between donor-hydrogen-acceptor in degrees."""
    if np:
        v1 = np.array(donor) - np.array(donor_h)
        v2 = np.array(acceptor) - np.array(donor_h)
        
        # Normalize vectors
        v1_norm = v1 / np.linalg.norm(v1)
        v2_norm = v2 / np.linalg.norm(v2)
        
        # Calculate angle
        cos_angle = np.clip(np.dot(v1_norm, v2_norm), -1.0, 1.0)
        return math.degrees(math.acos(cos_angle))
    else:
        # Pure Python implementation
        def subtract(a, b):
            return (a[0]-b[0], a[1]-b[1], a[2]-b[2])
        
        def normalize(v):
            length = math.sqrt(v[0]**2 + v[1]**2 + v[2]**2)
            return (v[0]/length, v[1]/length, v[2]/length)
        
        def dot(a, b):
            return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]
        
        v1 = normalize(subtract(donor, donor_h))
        v2 = normalize(subtract(acceptor, donor_h))
        
        cos_angle = max(-1.0, min(1.0, dot(v1, v2)))
        return math.degrees(math.acos(cos_angle))


def read_pdb_file(pdb_file: Path) -> List[Atom]:
    """Read PDB file and return list of atoms."""
    atoms = []
    with pdb_file.open('r', encoding='utf-8') as f:
        for line in f:
            if line.startswith(('ATOM', 'HETATM')):
                atoms.append(Atom(line))
    return atoms


def is_donor(atom: Atom) -> bool:
    """Check if atom can be a hydrogen bond donor."""
    # Common donor atoms: N, O with attached hydrogen
    if atom.element not in ['N', 'O']:
        return False
    
    # Check if atom name suggests it's attached to hydrogen
    donor_patterns = ['N', 'NH', 'NH1', 'NH2', 'NZ', 'NE', 'ND', 'OG', 'OH', 'SER', 'THR', 'TYR']
    for pattern in donor_patterns:
        if pattern in atom.atom_name:
            return True
    
    return False


def is_acceptor(atom: Atom) -> bool:
    """Check if atom can be a hydrogen bond acceptor."""
    # Common acceptor atoms: N, O, S
    if atom.element not in ['N', 'O', 'S']:
        return False
    
    # Check if atom name suggests it can accept hydrogen
    acceptor_patterns = ['N', 'O', 'OD1', 'OD2', 'OE1', 'OE2', 'OG', 'OH', 'SD', 'SG']
    for pattern in acceptor_patterns:
        if pattern in atom.atom_name:
            return True
    
    return False


def find_hydrogen_atoms(atoms: List[Atom], donor: Atom) -> List[Atom]:
    """Find hydrogen atoms attached to a donor atom."""
    hydrogens = []
    for atom in atoms:
        if atom.element == 'H' and atom.residue_id == donor.residue_id:
            # Check if hydrogen is close to donor (within 1.5 Å)
            if calculate_distance(donor.coord, atom.coord) < 1.5:
                hydrogens.append(atom)
    return hydrogens


def count_hydrogen_bonds(protein_atoms: List[Atom], ligand_atoms: List[Atom],
                        distance_cutoff: float = 3.5, angle_cutoff: float = 120.0) -> List[HydrogenBond]:
    """Count hydrogen bonds between protein and ligand."""
    hbonds = []
    
    # Create lookup dictionaries for faster access
    protein_by_id = {atom.atom_id: atom for atom in protein_atoms}
    ligand_by_id = {atom.atom_id: atom for atom in ligand_atoms}
    
    # Check protein donors to ligand acceptors
    for protein_atom in protein_atoms:
        if not is_donor(protein_atom):
            continue
            
        # Find hydrogens attached to this donor
        hydrogens = find_hydrogen_atoms(protein_atoms, protein_atom)
        
        for ligand_atom in ligand_atoms:
            if not is_acceptor(ligand_atom):
                continue
                
            for hydrogen in hydrogens:
                distance = calculate_distance(hydrogen.coord, ligand_atom.coord)
                if distance < distance_cutoff:
                    angle = calculate_angle(hydrogen.coord, protein_atom.coord, ligand_atom.coord)
                    if angle > angle_cutoff:
                        hbonds.append(HydrogenBond(protein_atom, ligand_atom, distance, angle))
    
    # Check ligand donors to protein acceptors
    for ligand_atom in ligand_atoms:
        if not is_donor(ligand_atom):
            continue
            
        # Find hydrogens attached to this donor
        hydrogens = find_hydrogen_atoms(ligand_atoms, ligand_atom)
        
        for protein_atom in protein_atoms:
            if not is_acceptor(protein_atom):
                continue
                
            for hydrogen in hydrogens:
                distance = calculate_distance(hydrogen.coord, protein_atom.coord)
                if distance < distance_cutoff:
                    angle = calculate_angle(hydrogen.coord, ligand_atom.coord, protein_atom.coord)
                    if angle > angle_cutoff:
                        hbonds.append(HydrogenBond(ligand_atom, protein_atom, distance, angle))
    
    return hbonds


def cluster_ligands(atoms: List[Atom], ligand_resname: str = "LIG", 
                    distance_cutoff: float = 2.0) -> Dict[int, List[Atom]]:
    """Cluster ligands based on spatial proximity."""
    ligand_residues = {}
    
    # Group atoms by residue
    for atom in atoms:
        if atom.residue_name == ligand_resname:
            if atom.residue_id not in ligand_residues:
                ligand_residues[atom.residue_id] = []
            ligand_residues[atom.residue_id].append(atom)
    
    # Simple clustering based on distance
    clusters = {}
    cluster_id = 0
    
    for res_id, ligand_atoms in ligand_residues.items():
        # Calculate center of mass
        com = [0.0, 0.0, 0.0]
        for atom in ligand_atoms:
            com[0] += atom.coord[0]
            com[1] += atom.coord[1]
            com[2] += atom.coord[2]
        com = [c/len(ligand_atoms) for c in com]
        
        # Check if this ligand belongs to an existing cluster
        assigned = False
        for existing_cluster_id, existing_cluster in clusters.items():
            existing_com = existing_cluster['com']
            distance = calculate_distance(com, existing_com)
            if distance < distance_cutoff:
                # Add to existing cluster
                existing_cluster['ligands'].append(res_id)
                # Update cluster center
                n = len(existing_cluster['ligands'])
                existing_com[0] = (existing_com[0] * (n-1) + com[0]) / n
                existing_com[1] = (existing_com[1] * (n-1) + com[1]) / n
                existing_com[2] = (existing_com[2] * (n-1) + com[2]) / n
                assigned = True
                break
        
        if not assigned:
            # Create new cluster
            clusters[cluster_id] = {
                'com': com,
                'ligands': [res_id]
            }
            cluster_id += 1
    
    return clusters


def calculate_interaction_energy(ligand_atoms: List[Atom], protein_atoms: List[Atom]) -> float:
    """Calculate simplified interaction energy between ligand and protein."""
    # This is a simplified calculation - in practice you'd use more sophisticated methods
    total_energy = 0.0
    
    for ligand_atom in ligand_atoms:
        for protein_atom in protein_atoms:
            distance = calculate_distance(ligand_atom.coord, protein_atom.coord)
            if distance < 0.1:  # Avoid division by zero
                continue
                
            # Simple Lennard-Jones-like potential
            # Using arbitrary parameters for demonstration
            sigma = 3.5  # Angstroms
            epsilon = 0.1  # kcal/mol
            
            r6 = (sigma / distance) ** 6
            r12 = r6 * r6
            
            # LJ potential: 4*epsilon*(r12 - r6)
            lj_energy = 4 * epsilon * (r12 - r6)
            
            # Simple electrostatic term (if both atoms have charges)
            # This is highly simplified - real calculations need proper charges
            if ligand_atom.element in ['N', 'O', 'S'] and protein_atom.element in ['N', 'O', 'S']:
                # Coulomb's law (simplified)
                q1 = 1.0 if ligand_atom.element == 'N' else -1.0
                q2 = 1.0 if protein_atom.element == 'N' else -1.0
                elec_energy = 332.0 * q1 * q2 / distance  # 332 is conversion factor
                total_energy += elec_energy * 0.1  # Scale down to reasonable range
            
            total_energy += lj_energy
    
    return total_energy


def calculate_wns_score(pdb_file: Path, trajectory_file: Optional[Path] = None,
                        ligand_resname: str = "LIG") -> List[Dict]:
    """
    Calculate WnS score for ligands.
    WnS Score = interaction energy + hydrogen bond contribution
    """
    print(f"Analyzing hydrogen bonds in {pdb_file}")
    
    # Read structure
    atoms = read_pdb_file(pdb_file)
    
    # Separate protein and ligand atoms
    protein_atoms = [atom for atom in atoms if atom.residue_name != ligand_resname]
    ligand_residues = {}
    for atom in atoms:
        if atom.residue_name == ligand_resname:
            if atom.residue_id not in ligand_residues:
                ligand_residues[atom.residue_id] = []
            ligand_residues[atom.residue_id].append(atom)
    
    if not ligand_residues:
        print(f"No ligands with residue name '{ligand_resname}' found")
        return []
    
    print(f"Found {len(ligand_residues)} ligand residues")
    
    # Cluster ligands
    clusters = cluster_ligands(atoms, ligand_resname)
    print(f"Found {len(clusters)} ligand clusters")
    
    scores = []
    
    for cluster_id, cluster_info in clusters.items():
        # Get representative ligand (closest to cluster center)
        cluster_com = cluster_info['com']
        best_ligand_id = None
        best_distance = float('inf')
        
        for ligand_id in cluster_info['ligands']:
            ligand_atoms = ligand_residues[ligand_id]
            ligand_com = [0.0, 0.0, 0.0]
            for atom in ligand_atoms:
                ligand_com[0] += atom.coord[0]
                ligand_com[1] += atom.coord[1]
                ligand_com[2] += atom.coord[2]
            ligand_com = [c/len(ligand_atoms) for c in ligand_com]
            
            distance = calculate_distance(ligand_com, cluster_com)
            if distance < best_distance:
                best_distance = distance
                best_ligand_id = ligand_id
        
        if best_ligand_id is None:
            continue
            
        ligand_atoms = ligand_residues[best_ligand_id]
        
        # Calculate interaction energy
        e_inter = calculate_interaction_energy(ligand_atoms, protein_atoms)
        
        # Count hydrogen bonds
        hbonds = count_hydrogen_bonds(protein_atoms, ligand_atoms)
        n_hbonds = len(hbonds)
        
        # Only process ligands that have hydrogen bonds
        if n_hbonds == 0:
            print(f"Ligand {best_ligand_id}: No hydrogen bonds found, skipping")
            continue
        
        # Calculate WnS score (lower is better)
        # Energy term (negative for favorable) + hydrogen bond bonus
        final_score = e_inter - (n_hbonds * 2.0)
        
        # Prepare hydrogen bond details
        hbond_details = []
        for hbond in hbonds:
            hbond_details.append({
                'donor_residue': f"{hbond.donor.residue_name}{hbond.donor.residue_id}",
                'donor_atom': hbond.donor.atom_name,
                'acceptor_residue': f"{hbond.acceptor.residue_name}{hbond.acceptor.residue_id}",
                'acceptor_atom': hbond.acceptor.atom_name,
                'distance': hbond.distance,
                'angle': hbond.angle
            })
        
        scores.append({
            'ligand_id': best_ligand_id,
            'cluster_id': cluster_id,
            'score': final_score,
            'e_inter': e_inter,
            'n_hbonds': n_hbonds,
            'hbond_details': hbond_details,
            'coords': cluster_com
        })
        
        print(f"Ligand {best_ligand_id}: Score={final_score:.2f}, H-bonds={n_hbonds}, E_inter={e_inter:.2f}")
    
    # Sort by score (lower is better)
    return sorted(scores, key=lambda x: x['score'])





def save_results_to_csv(scores: List[Dict], output_csv: str = "hbond_analysis.csv") -> None:
    """Save analysis results to CSV file."""
    print(f"Saving results to {output_csv}")
    
    with open(output_csv, 'w', newline='', encoding='utf-8') as csvfile:
        fieldnames = ['ligand_id', 'cluster_id', 'score', 'e_inter', 'n_hbonds', 
                     'hbond_details', 'coords_x', 'coords_y', 'coords_z']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        
        writer.writeheader()
        for item in scores:
            # Convert hbond_details to string for CSV
            hbond_str = '; '.join([
                f"{detail['donor_residue']}({detail['donor_atom']})-{detail['acceptor_residue']}({detail['acceptor_atom']}) "
                f"d={detail['distance']:.2f}Å a={detail['angle']:.1f}°"
                for detail in item['hbond_details']
            ])
            
            writer.writerow({
                'ligand_id': item['ligand_id'],
                'cluster_id': item['cluster_id'],
                'score': item['score'],
                'e_inter': item['e_inter'],
                'n_hbonds': item['n_hbonds'],
                'hbond_details': hbond_str,
                'coords_x': item['coords'][0],
                'coords_y': item['coords'][1],
                'coords_z': item['coords'][2]
            })
    
    print(f"Results saved to {output_csv}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("pdb_file", type=Path, help="PDB file with protein-ligand complex")
    parser.add_argument("-t", "--trajectory", type=Path, help="Optional trajectory file for analysis")
    parser.add_argument("-l", "--ligand", default="LIG", help="Ligand residue name (default: LIG)")
    parser.add_argument("-o", "--output", default="view_hbonds", help="Output file prefix")
    parser.add_argument("-d", "--distance", type=float, default=3.5, 
                       help="Hydrogen bond distance cutoff (Å)")
    parser.add_argument("-a", "--angle", type=float, default=120.0,
                       help="Hydrogen bond angle cutoff (degrees)")
    parser.add_argument("--csv", action="store_true", help="Save results to CSV file")
    
    args = parser.parse_args()
    
    if not args.pdb_file.exists():
        raise FileNotFoundError(f"PDB file not found: {args.pdb_file}")
    
    # Calculate WnS scores
    scores = calculate_wns_score(args.pdb_file, args.trajectory, args.ligand)
    
    if not scores:
        print("No ligands found for analysis")
        return
    
    # Save to CSV if requested
    if args.csv:
        csv_file = f"{args.output}.csv"
        save_results_to_csv(scores, csv_file)
    
    # Print summary
    print("\n=== Hydrogen Bond Analysis Summary ===")
    print(f"Total ligands with hydrogen bonds: {len(scores)}")
    print(f"Total hydrogen bonds found: {sum(item['n_hbonds'] for item in scores)}")
    
    print("\nTop ligands by WnS score:")
    for i, item in enumerate(scores):
        print(f"{i+1}. Ligand {item['ligand_id']}: Score={item['score']:.2f}, "
              f"H-bonds={item['n_hbonds']}, Energy={item['e_inter']:.2f}")
    
    print(f"\nResults saved to: {args.output}.csv")
    print("Only ligands with hydrogen bonds are included in the output")


if __name__ == "__main__":
    main()