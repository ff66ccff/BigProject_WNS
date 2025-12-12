#!/usr/bin/env python3
"""Generate ligand topology using OpenBabel and ParmEd."""

import os
import sys
from pathlib import Path

def generate_ligand_itp(mol2_file, output_itp):
    """Generate a simple ITP file from MOL2 using OpenBabel and ParmEd."""
    
    # Convert MOL2 to PDB with proper atom names
    pdb_file = mol2_file.with_suffix('.pdb')
    os.system(f"obabel -imol2 {mol2_file} -opdb -O {pdb_file} -h")
    
    # Read the PDB file to get atom information
    atoms = []
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                atom_num = int(line[6:11])
                atom_name = line[12:16].strip()
                res_name = line[17:20].strip()
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                element = line[76:78].strip() if len(line) > 76 else atom_name[0]
                atoms.append({
                    'num': atom_num,
                    'name': atom_name,
                    'res_name': res_name,
                    'element': element,
                    'x': x, 'y': y, 'z': z
                })
    
    # Generate bonds based on distance (simplified approach)
    bonds = []
    for i in range(len(atoms)):
        for j in range(i+1, len(atoms)):
            dx = atoms[i]['x'] - atoms[j]['x']
            dy = atoms[i]['y'] - atoms[j]['y']
            dz = atoms[i]['z'] - atoms[j]['z']
            dist = (dx*dx + dy*dy + dz*dz) ** 0.5
            
            # Simple bond detection based on distance and element types
            if dist < 1.6:  # Typical bond length threshold
                bonds.append((i+1, j+1))  # GROMACS uses 1-based indexing
    
    # Write ITP file
    with open(output_itp, 'w') as f:
        f.write("[ moleculetype ]\n")
        f.write("; Name            nrexcl\n")
        f.write("UNL              3\n\n")
        
        f.write("[ atoms ]\n")
        f.write("; nr  type  resnr  residu atom  cgnr  charge   mass\n")
        
        # Simple atom type mapping
        type_map = {'H': 'H', 'C': 'C', 'N': 'N', 'O': 'O', 'S': 'S'}
        
        for i, atom in enumerate(atoms):
            atom_type = type_map.get(atom['element'], 'C')
            # Simple charge assignment (neutral molecule)
            charge = 0.0
            if atom['element'] == 'H':
                charge = 0.0
            elif atom['element'] == 'C':
                charge = -0.1
            elif atom['element'] == 'N':
                charge = 0.1
            elif atom['element'] == 'O':
                charge = -0.2
                
            mass = {'H': 1.008, 'C': 12.011, 'N': 14.007, 'O': 15.999, 'S': 32.066}
            
            f.write(f"{i+1:4d} {atom_type:5s}     1      UNL    {atom['name']:4s}    {i+1:4d}    {charge:8.4f}   {mass.get(atom['element'], 12.011):8.5f}\n")
        
        f.write("\n[ bonds ]\n")
        f.write("; ai   aj   funct   length   force.c.\n")
        for bond in bonds:
            f.write(f"{bond[0]:4d} {bond[1]:4d}    1\n")
        
        f.write("\n[ system ]\n")
        f.write("; Name\n")
        f.write("Ligand\n\n")
        
        f.write("[ molecules ]\n")
        f.write("; Compound        #mols\n")
        f.write("UNL              1\n")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python generate_ligand_top.py input.mol2 output.itp")
        sys.exit(1)
    
    mol2_file = Path(sys.argv[1])
    output_itp = Path(sys.argv[2])
    
    generate_ligand_itp(mol2_file, output_itp)
    print(f"Generated {output_itp}")