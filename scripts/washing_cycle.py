#!/usr/bin/env python3
"""Washing cycle for Wrap 'n' Shake: MD with ligand displacement analysis."""

from __future__ import annotations

import argparse
import math
import os
import re
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Tuple, Set


class AtomRecord:
    """Represents an atom in a GRO file."""
    
    def __init__(self, line: str) -> None:
        self.line = line.rstrip('\n')
        if len(line.split()) >= 6:
            self.atom_id = int(line.split()[0])
            self.atom_name = line.split()[1]
            self.residue_name = line.split()[2]
            self.residue_id = int(line.split()[3])
            self.coord = (
                float(line.split()[4]),
                float(line.split()[5]),
                float(line.split()[6]),
            )
        else:
            self.atom_id = 0
            self.atom_name = ""
            self.residue_name = ""
            self.residue_id = 0
            self.coord = (0.0, 0.0, 0.0)


def calculate_distance(coord1: Tuple[float, float, float], 
                      coord2: Tuple[float, float, float]) -> float:
    """Calculate Euclidean distance between two coordinates."""
    return math.sqrt(
        (coord1[0] - coord2[0])**2 + 
        (coord1[1] - coord2[1])**2 + 
        (coord1[2] - coord2[2])**2
    )


def calculate_rmsd(coords1: List[Tuple[float, float, float]], 
                  coords2: List[Tuple[float, float, float]]) -> float:
    """Calculate RMSD between two sets of coordinates."""
    if len(coords1) != len(coords2):
        return float('inf')
    
    sum_sq_diff = 0.0
    for c1, c2 in zip(coords1, coords2):
        sum_sq_diff += (c1[0] - c2[0])**2 + (c1[1] - c2[1])**2 + (c1[2] - c2[2])**2
    
    return math.sqrt(sum_sq_diff / len(coords1))


def read_gro_file(gro_file: Path) -> Tuple[List[AtomRecord], List[str]]:
    """Read GRO file and return atom records and header/footer."""
    lines = gro_file.read_text(encoding='utf-8').splitlines()
    
    if len(lines) < 2:
        raise ValueError(f"Invalid GRO file: {gro_file}")
    
    header = lines[0]
    atom_count = int(lines[1])
    atom_lines = lines[2:2+atom_count]
    footer = lines[2+atom_count] if len(lines) > 2+atom_count else ""
    
    atoms = [AtomRecord(line) for line in atom_lines]
    return atoms, [header, str(atom_count), footer]


def write_gro_file(gro_file: Path, atoms: List[AtomRecord], 
                  header_footer: List[str]) -> None:
    """Write atoms to GRO file with given header and footer."""
    with gro_file.open('w', encoding='utf-8') as f:
        f.write(header_footer[0] + '\n')
        f.write(str(len(atoms)) + '\n')
        for atom in atoms:
            f.write(atom.line + '\n')
        f.write(header_footer[2] + '\n')


def get_ligand_residues(atoms: List[AtomRecord], ligand_resname: str) -> Dict[int, List[AtomRecord]]:
    """Group atoms by ligand residue ID."""
    ligand_residues = {}
    
    for atom in atoms:
        if atom.residue_name == ligand_resname:
            if atom.residue_id not in ligand_residues:
                ligand_residues[atom.residue_id] = []
            ligand_residues[atom.residue_id].append(atom)
    
    return ligand_residues


def calculate_ligand_displacements(initial_atoms: List[AtomRecord], 
                                 final_atoms: List[AtomRecord],
                                 ligand_resname: str) -> Dict[int, float]:
    """Calculate displacement for each ligand residue."""
    initial_ligands = get_ligand_residues(initial_atoms, ligand_resname)
    final_ligands = get_ligand_residues(final_atoms, ligand_resname)
    
    displacements = {}
    
    for res_id, initial_res_atoms in initial_ligands.items():
        if res_id not in final_ligands:
            continue  # Ligand disappeared
        
        final_res_atoms = final_ligands[res_id]
        
        # Calculate center of mass displacement
        initial_com = [0.0, 0.0, 0.0]
        final_com = [0.0, 0.0, 0.0]
        
        for atom in initial_res_atoms:
            initial_com[0] += atom.coord[0]
            initial_com[1] += atom.coord[1]
            initial_com[2] += atom.coord[2]
        
        for atom in final_res_atoms:
            final_com[0] += atom.coord[0]
            final_com[1] += atom.coord[1]
            final_com[2] += atom.coord[2]
        
        n_atoms = len(initial_res_atoms)
        initial_com = [c/n_atoms for c in initial_com]
        final_com = [c/n_atoms for c in final_com]
        
        displacement = calculate_distance(tuple(initial_com), tuple(final_com))
        displacements[res_id] = displacement
    
    return displacements


def update_topology_file(topol_file: Path, ligand_resname: str, 
                        removed_residues: Set[int]) -> None:
    """Update topology file to remove washed-away ligands."""
    lines = topol_file.read_text(encoding='utf-8').splitlines()
    
    # Find molecules section
    molecules_start = -1
    for i, line in enumerate(lines):
        if '[ molecules ]' in line:
            molecules_start = i + 2  # Skip header and blank line
            break
    
    if molecules_start == -1:
        raise ValueError(f"Could not find [ molecules ] section in {topol_file}")
    
    # Update ligand count
    ligand_count_updated = False
    for i in range(molecules_start, len(lines)):
        if ligand_resname in lines[i] and not lines[i].strip().startswith(';'):
            parts = lines[i].split()
            if len(parts) >= 2:
                try:
                    original_count = int(parts[1])
                    new_count = max(0, original_count - len(removed_residues))
                    parts[1] = str(new_count)
                    lines[i] = ' '.join(parts)
                    ligand_count_updated = True
                    print(f"Updated ligand count: {original_count} -> {new_count}")
                    break
                except ValueError:
                    continue
    
    if not ligand_count_updated:
        raise ValueError(f"Could not find ligand {ligand_resname} in molecules section")
    
    # Write updated topology
    topol_file.write_text('\n'.join(lines) + '\n', encoding='utf-8')
    print(f"Updated {topol_file}: removed {len(removed_residues)} {ligand_resname} molecules")


def run_gromacs_command(gmx_exe: str, args: List[str], cwd: Path, 
                       input_text: str = "") -> str:
    """Run GROMACS command and return output."""
    cmd = [gmx_exe] + args
    print(f"Running: {' '.join(cmd)}")
    
    result = subprocess.run(
        cmd, 
        cwd=cwd, 
        input=input_text, 
        text=True, 
        capture_output=True
    )
    
    if result.returncode != 0:
        raise RuntimeError(f"GROMACS command failed: {result.stderr}")
    
    return result.stdout


def regenerate_index_files(work_dir: Path, gmx_exe: str, gro_file: str = "npt.gro") -> None:
    """Regenerate GROMACS index files to prevent index mismatch."""
    print("Regenerating index files...")
    
    # Create new index file
    run_gromacs_command(
        gmx_exe,
        ["make_ndx", "-f", gro_file, "-o", "index.ndx"],
        work_dir,
        input_text="q\n"  # Quit after creating default groups
    )
    
    print("Index files regenerated")


def validate_annealing_mdp(mdp_file: Path, cycle_time: float) -> None:
    """Validate that annealing.mdp parameters are consistent with cycle time."""
    if not mdp_file.exists():
        raise FileNotFoundError(f"Annealing MDP file not found: {mdp_file}")
    
    content = mdp_file.read_text(encoding='utf-8')
    
    # Extract nsteps and dt
    nsteps_match = re.search(r'nsteps\s*=\s*(\d+)', content)
    dt_match = re.search(r'dt\s*=\s*([\d.]+)', content)
    
    if not nsteps_match or not dt_match:
        raise ValueError("Could not find nsteps or dt in annealing.mdp")
    
    nsteps = int(nsteps_match.group(1))
    dt = float(dt_match.group(1))
    
    # Calculate total simulation time in ps
    total_time_ps = nsteps * dt * 1000  # dt is in ps, so multiply by 1000 for ps
    
    # Extract annealing time points
    annealing_time_match = re.search(r'annealing-time\s*=\s*([\d\s.]+)', content)
    if not annealing_time_match:
        raise ValueError("Could not find annealing-time in annealing.mdp")
    
    annealing_times = [float(t) for t in annealing_time_match.group(1).split()]
    max_annealing_time = max(annealing_times)
    
    # Check if annealing time covers the full simulation
    if max_annealing_time > total_time_ps:
        raise ValueError(
            f"Annealing time ({max_annealing_time} ps) exceeds total simulation time "
            f"({total_time_ps} ps). Please adjust nsteps or dt."
        )
    
    # Check if cycle time matches
    expected_time_ns = total_time_ps / 1000  # Convert to ns
    if abs(expected_time_ns - cycle_time) > 0.01:
        print(f"WARNING: Expected cycle time ({cycle_time} ns) differs from MDP time ({expected_time_ns} ns)")
    
    print(f"Annealing MDP validation passed: {total_time_ps} ps total time")


def washing_cycle(work_dir: Path, gmx_exe: str = "gmx", 
                ligand_resname: str = "LIG", 
                displacement_cutoff: float = 6.0,
                cycle_time: float = 1.0) -> None:
    """Run washing cycle with MD and ligand displacement analysis.
    
    Args:
        work_dir: Working directory with GROMACS files
        gmx_exe: GROMACS executable name
        ligand_resname: Residue name for ligands
        displacement_cutoff: Displacement threshold in Angstroms
        cycle_time: MD simulation time per cycle in nanoseconds
    """
    print(f"Starting washing cycle in {work_dir}")
    
    # Check required files
    required_files = ["topol.top", "npt.gro", "annealing.mdp"]
    for filename in required_files:
        if not (work_dir / filename).exists():
            raise FileNotFoundError(f"Required file not found: {work_dir / filename}")
    
    # Validate annealing MDP parameters
    validate_annealing_mdp(work_dir / "annealing.mdp", cycle_time)
    
    # Regenerate index files to prevent index mismatch
    regenerate_index_files(work_dir, gmx_exe)
    
    # Read initial structure
    initial_atoms, header_footer = read_gro_file(work_dir / "npt.gro")
    initial_ligands = get_ligand_residues(initial_atoms, ligand_resname)
    
    if not initial_ligands:
        print("No ligands found in the system")
        return
    
    print(f"Found {len(initial_ligands)} ligand residues")
    
    # Run annealing MD
    print(f"Running {cycle_time} ns annealing simulation...")
    run_gromacs_command(
        gmx_exe, 
        ["grompp", "-f", "annealing.mdp", "-c", "npt.gro", "-p", "topol.top", "-o", "anneal.tpr", "-maxwarn", "1"],
        work_dir
    )
    
    run_gromacs_command(
        gmx_exe, 
        ["mdrun", "-deffnm", "anneal"],
        work_dir
    )
    
    # Extract final frame - calculate correct time point
    final_time_ps = int(cycle_time * 1000)  # Convert ns to ps
    print(f"Extracting final frame at {final_time_ps} ps...")
    run_gromacs_command(
        gmx_exe,
        ["trjconv", "-s", "anneal.tpr", "-f", "anneal.xtc", "-o", "final.gro", "-dump", str(final_time_ps)],
        work_dir,
        input_text="0\n"  # Select system
    )
    
    # Read final structure
    final_atoms, _ = read_gro_file(work_dir / "final.gro")
    
    # Calculate displacements
    displacements = calculate_ligand_displacements(initial_atoms, final_atoms, ligand_resname)
    
    # Identify washed-away ligands
    washed_residues = {res_id for res_id, disp in displacements.items() if disp > displacement_cutoff}
    
    print(f"Ligand displacements:")
    for res_id, disp in displacements.items():
        status = "WASHED AWAY" if disp > displacement_cutoff else "OK"
        print(f"  Residue {res_id}: {disp:.2f} Å [{status}]")
    
    if not washed_residues:
        print("No ligands washed away in this cycle")
        return
    
    print(f"Removing {len(washed_residues)} washed-away ligands...")
    
    # Filter atoms to remove washed ligands
    filtered_atoms = [
        atom for atom in final_atoms 
        if not (atom.residue_name == ligand_resname and atom.residue_id in washed_residues)
    ]
    
    # Write filtered structure
    write_gro_file(work_dir / "filtered.gro", filtered_atoms, header_footer)
    
    # Update topology
    update_topology_file(work_dir / "topol.top", ligand_resname, washed_residues)
    
    # Update atom numbers in filtered structure
    for i, atom in enumerate(filtered_atoms):
        atom.atom_id = i + 1
        # Update the line with new atom ID
        parts = atom.line.split()
        parts[0] = str(i + 1)
        atom.line = ' '.join(parts)
    
    # Update header with new atom count
    header_footer[1] = str(len(filtered_atoms))
    
    write_gro_file(work_dir / "npt.gro", filtered_atoms, header_footer)
    
    # CRITICAL: Regenerate index files after topology change
    regenerate_index_files(work_dir, gmx_exe)
    
    print(f"Washing cycle completed. {len(washed_residues)} ligands removed.")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("work_dir", type=Path, help="Working directory with GROMACS files")
    parser.add_argument("-g", "--gmx", default="gmx", help="GROMACS executable")
    parser.add_argument("-l", "--ligand", default="LIG", help="Ligand residue name")
    parser.add_argument("-c", "--cutoff", type=float, default=6.0, 
                       help="Displacement cutoff in Angstroms")
    parser.add_argument("-t", "--time", type=float, default=1.0,
                       help="MD simulation time per cycle in nanoseconds")
    
    args = parser.parse_args()
    
    if not args.work_dir.exists():
        raise FileNotFoundError(f"Working directory not found: {args.work_dir}")
    
    washing_cycle(
        args.work_dir, 
        args.gmx, 
        args.ligand, 
        args.cutoff,
        args.time
    )


if __name__ == "__main__":
    main()