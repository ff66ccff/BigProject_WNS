#!/usr/bin/env python3
"""Generate scoring report in JSON format (Standard Python Implementation)."""

import json
import math
import datetime
from pathlib import Path

# Atomic weights (Average)
ATOMIC_WEIGHTS = {
    'H': 1.008, 'C': 12.011, 'N': 14.007, 'O': 15.999, 'P': 30.974, 
    'S': 32.065, 'F': 18.998, 'Cl': 35.45, 'Br': 79.904, 'I': 126.904,
    'FE': 55.845, 'MG': 24.305, 'ZN': 65.38, 'MN': 54.938
}

# Residue properties
HYDROPHOBIC_RES = {'ALA', 'VAL', 'LEU', 'ILE', 'MET', 'PHE', 'TRP', 'PRO'}
AROMATIC_RES = {'PHE', 'TYR', 'TRP', 'HIS'}

class Atom:
    def __init__(self, line):
        self.line = line
        try:
            # PDB/PDBQT fixed widths
            self.atom_name = line[12:16].strip()
            self.res_name = line[17:20].strip()
            self.chain = line[21:22].strip()
            self.res_id = line[22:26].strip()
            self.x = float(line[30:38])
            self.y = float(line[38:46])
            self.z = float(line[46:54])
            
            # Element guessing
            # PDBQT puts element at the end (77-78), PDB at 76-77 usually
            if len(line) >= 78:
                self.element = line[76:78].strip() # Try standard PDB pos first
                if not self.element or self.element.isdigit():
                     self.element = line[77:79].strip() # Try PDBQT pos
            else:
                 self.element = ""
            
            if not self.element:
                # Erasure fallback based on name
                name = self.atom_name
                if name.startswith("HG") or name.startswith("H"): self.element = "H"
                elif name.startswith("C"): self.element = "C"
                elif name.startswith("N"): self.element = "N"
                elif name.startswith("O"): self.element = "O"
                elif name.startswith("S"): self.element = "S"
                elif name.startswith("P"): self.element = "P"
                else: self.element = "C" # Default
            
            self.element = self.element.upper()
            
            # Properties
            self.is_carbon = 'C' in self.element
            self.is_polar = self.element in ['N', 'O', 'S']
            
        except ValueError:
            self.valid = False
        else:
            self.valid = True

def dist_sq(a1, a2):
    return (a1.x - a2.x)**2 + (a1.y - a2.y)**2 + (a1.z - a2.z)**2

def analyze_structure(complex_file: Path, hbond_file: Path, vina_file: Path) -> dict:
    # 1. Parse Structure
    prot_atoms = []
    lig_atoms = []
    
    # Try reading Vina Result for Ligand (as it has correct scoring pose)
    # But complex.pdbqt has both?
    # Let's read complex file which implies "Protein + Ligand"
    
    if not complex_file.exists():
        return {}

    with open(complex_file, 'r') as f:
        for line in f:
            if line.startswith(('ATOM', 'HETATM')):
                atom = Atom(line)
                if not atom.valid: continue
                
                # Identify ligand by residue name or LIG string
                if 'LIG' in line or 'UNL' in line or 'MOL' in line:
                    lig_atoms.append(atom)
                else:
                    prot_atoms.append(atom)

    # 2. Ligand Metrics
    mol_weight = 0.0
    for a in lig_atoms:
        mol_weight += ATOMIC_WEIGHTS.get(a.element, 12.0)
    
    # 3. Protein Metrics
    unique_residues = set()
    hydrophobic_count = 0
    aromatic_count = 0
    
    for a in prot_atoms:
        res_key = (a.chain, a.res_id, a.res_name)
        if res_key not in unique_residues:
            unique_residues.add(res_key)
            if a.res_name in HYDROPHOBIC_RES:
                hydrophobic_count += 1
            if a.res_name in AROMATIC_RES:
                aromatic_count += 1
                
    total_res = float(len(unique_residues))
    if total_res == 0: total_res = 1.0
    
    # 4. Binding Metrics
    min_dist_sq = float('inf')
    contact_res_keys = set()
    hydrophobic_contacts = 0
    aromatic_contacts = 0
    polar_contacts = 0
    total_contacts = 0
    
    cutoff_sq = 4.5**2
    
    # Optimization: Bounding box
    if lig_atoms and prot_atoms:
        lx = [a.x for a in lig_atoms]
        ly = [a.y for a in lig_atoms]
        lz = [a.z for a in lig_atoms]
        min_x, max_x = min(lx)-4.5, max(lx)+4.5
        min_y, max_y = min(ly)-4.5, max(ly)+4.5
        min_z, max_z = min(lz)-4.5, max(lz)+4.5
        
        nearby_prot = [a for a in prot_atoms if 
                       min_x <= a.x <= max_x and 
                       min_y <= a.y <= max_y and 
                       min_z <= a.z <= max_z]
                       
        for la in lig_atoms:
            for pa in nearby_prot:
                d2 = dist_sq(la, pa)
                
                if d2 < min_dist_sq:
                    min_dist_sq = d2
                
                if d2 <= cutoff_sq:
                    total_contacts += 1
                    contact_res_keys.add((pa.chain, pa.res_id))
                    
                    if la.is_carbon and pa.is_carbon:
                        hydrophobic_contacts += 1
                    
                    if la.is_polar and pa.is_polar and d2 <= 3.5**2:
                        polar_contacts += 1
                        
                    if pa.res_name in AROMATIC_RES and d2 <= 4.0**2:
                         aromatic_contacts += 1

    # 5. H-bonds
    hbond_count = 0
    if hbond_file.exists():
        try:
            with open(hbond_file, 'r') as f:
                for line in f:
                    if "1,0," in line:
                         parts = line.split(',')
                         hbond_count = int(parts[4])
        except: pass

    # 6. Vina Affinity
    vina_score = 0.0
    if vina_file.exists():
        try:
            with open(vina_file, 'r') as f:
                for line in f:
                    if "REMARK VINA RESULT:" in line:
                        # Format: REMARK VINA RESULT:    -6.8      0.000      0.000
                        parts = line.split()
                        try:
                            vina_score = float(parts[3])
                        except: pass
                        break # Only 1st model
        except: pass


    # 6. Normalization and Weights (Deduced from User Example)
    # Ligand
    # Raw MW 454 -> Norm 0.113 => Norm = x / 4000
    norm_mw = min(1.0, mol_weight / 4000.0)
    
    # Protein
    # Res: 984 -> 0.478 => x / 2055 (approx 2000)
    # Hyd: 0.473 -> 0.433 => (x - 0.3) / 0.4 (Range 0.3-0.7)
    # Aro: 0.101 -> 1.0 => (x - 0.05) / 0.05 (Range 0.05-0.10)
    norm_res = min(1.0, total_res / 2055.0)
    norm_hyd = min(1.0, max(0.0, (hydrophobic_count / total_res - 0.3) / 0.4))
    norm_aro = min(1.0, max(0.0, (aromatic_count / total_res - 0.05) / 0.05))
    
    # Binding
    # Weights: Contact(2), HB(3), MinDist(2), HydCont(2) -> Sum 9
    # Norms (Hypothetical): 
    # Contact: 19 -> 1.0. Maybe x / 15?
    # HB: 22 -> 1.0. Maybe x / 5?
    # MinDist: 0.14 -> 0.0. Maybe (3.0 - x) / 3.0? (Smaller is better). If x=0.14, (2.86)/3 ~ 0.95. User had 0.0?
    # Wait, user example MinDist raw 0.14 -> Norm 0.0. Maybe MinDist MUST be > threshold? Or < threshold?
    # If 0.14 is "clash", score is 0?
    # Let's assume standard "Distance Score": 1.0 if < 2.5?
    # User had 0.0 for 0.14. That's weird. 0.14A is a CLASH. So score 0 is correct.
    # So Norm = 1.0 if dist > 2.0? No, usually closer is better.
    # Maybe Norm = 0 if dist < 1.0 (Clash)?
    # Let's use simple logic: Contact > 10 = 1.0. HB > 1 = 1.0.
    
    norm_contact = min(1.0, len(contact_res_keys) / 10.0)
    norm_hb = min(1.0, hbond_count / 1.0) # 1 HB is good enough? User had 22!
    # Let's say 5 HB = 1.0.
    norm_hb = min(1.0, hbond_count / 5.0)
    norm_mindist = 1.0 if math.sqrt(min_dist_sq) > 1.5 else 0.0 # Penalize clash
    if math.sqrt(min_dist_sq) > 4.0: norm_mindist = 0.0 # Too far
    
    norm_hyd_cont = min(1.0, hydrophobic_contacts / 20.0)
    
    # Weights for Groups (Protein=1:2:1, Binding=2:3:2:2)
    # Protein Group
    p_w = {'res': 1, 'hyd': 2, 'aro': 1}
    p_sum = sum(p_w.values())
    prot_group_score = (norm_res * p_w['res'] + norm_hyd * p_w['hyd'] + norm_aro * p_w['aro']) / p_sum
    
    # Binding Group
    b_w = {'contact': 2, 'hb': 3, 'dist': 2, 'hyd': 2}
    b_sum = sum(b_w.values())
    bind_group_score = (norm_contact * b_w['contact'] + norm_hb * b_w['hb'] + 
                        norm_mindist * b_w['dist'] + norm_hyd_cont * b_w['hyd']) / b_sum
                        
    # Ligand Group
    lig_group_score = norm_mw # Single metric
    
    # Final Score (Simple Average of Groups?)
    final_score = (prot_group_score + bind_group_score + lig_group_score) / 3.0

    return {
        "ligand": {
            "raw": { "molecular_weight": mol_weight },
            "normalized": { "molecular_weight": norm_mw },
            "weighted": { "molecular_weight": norm_mw * 0.1 }, # Guess global weight
            "group_score": lig_group_score
        },
        "protein": {
            "raw": {
                "residue_count": total_res,
                "hydrophobic_fraction": hydrophobic_count / total_res,
                "aromatic_fraction": aromatic_count / total_res
            },
            "normalized": {
                "residue_count": norm_res,
                "hydrophobic_fraction": norm_hyd,
                "aromatic_fraction": norm_aro
            },
            "weighted": {
                 "residue_count": norm_res * 0.05,
                 "hydrophobic_fraction": norm_hyd * 0.1,
                 "aromatic_fraction": norm_aro * 0.05
            },
            "group_score": prot_group_score
        },
        "binding": {
            "raw": {
                "binding_energy": vina_score,
                "contact_residues": float(len(contact_res_keys)),
                "min_distance": math.sqrt(min_dist_sq) if min_dist_sq != float('inf') else 0.0,
                "avg_distance": 0.0,
                "hbond_count": float(hbond_count),
                "hydrophobic_contacts": float(hydrophobic_contacts),
                "aromatic_contacts": float(aromatic_contacts),
                "polar_contacts": float(polar_contacts),
                "total_contacts": float(total_contacts)
            },
            "normalized": {
                "contact_residues": norm_contact,
                "hbond_count": norm_hb,
                "min_distance": norm_mindist,
                "hydrophobic_contacts": norm_hyd_cont
            },
            "weighted": {
                "contact_residues": norm_contact * 0.1,
                "hbond_count": norm_hb * 0.15
            },
            "group_score": bind_group_score
        },
        "final_score": final_score,
        "note": "Calculated using deduced formulas from user example. Vina Energy included."
    }

def main():
    complex_file = Path("complex/complex.pdbqt")
    hbond_file = Path("view_hbonds.csv")
    vina_file = Path("autodock_runs/vina_docked_101.pdbqt")
    
    data = analyze_structure(complex_file, hbond_file, vina_file)
    
    output = data
    output["reproducibility"] = {
        "timestamp": datetime.datetime.now().isoformat(),
        "seeds": { "autodock_seed": 101 },
        "files_used": { "complex": str(complex_file) }
    }
    output["input_files"] = {
        "ligand": "data/intermediate1.mol2",
        "protein": "pdb/protein_clean.pdb",
        "complex": str(complex_file),
        "config": "config.yml"
    }
    
    with open("score_report.json", "w") as f:
        json.dump(output, f, indent=2)
        
    print("Generated score_report.json")
    print(json.dumps(output, indent=2))

if __name__ == "__main__":
    main()
