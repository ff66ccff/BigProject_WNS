#!/usr/bin/env python3
"""Extract best pose from Vina output and calculate hydrogen bonds."""

from pathlib import Path
import glob
import re

# Extract the best pose from each Vina output
def extract_best_poses(output_dir: Path):
    """Extract best pose (MODEL 1) from each Vina docking result."""
    results = []
    
    for pdbqt_file in sorted(output_dir.glob("vina_docked_*.pdbqt")):
        with open(pdbqt_file, 'r') as f:
            content = f.read()
        
        # Extract first MODEL
        match = re.search(r'MODEL 1\n(.*?)ENDMDL', content, re.DOTALL)
        if match:
            # Get affinity from REMARK
            affinity_match = re.search(r'VINA RESULT:\s+([-\d.]+)', match.group(0))
            affinity = float(affinity_match.group(1)) if affinity_match else 0.0
            
            results.append({
                'file': pdbqt_file.name,
                'affinity': affinity,
                'content': match.group(0)
            })
    
    # Sort by affinity (lower is better)
    results.sort(key=lambda x: x['affinity'])
    return results


def pdbqt_to_pdb(pdbqt_content: str) -> str:
    """Convert PDBQT content to PDB format."""
    pdb_lines = []
    for line in pdbqt_content.split('\n'):
        if line.startswith(('ATOM', 'HETATM')):
            # Convert PDBQT to PDB (remove charge and type columns)
            pdb_line = line[:66]
            # Replace UNL with LIG for ligand residue name
            pdb_line = pdb_line[:17] + 'LIG' + pdb_line[20:]
            pdb_lines.append(pdb_line)
    return '\n'.join(pdb_lines)


def main():
    output_dir = Path("autodock_runs")
    complex_dir = Path("complex")
    complex_dir.mkdir(exist_ok=True)
    
    # Extract best poses
    print("=== Extracting Best Poses from Vina Results ===")
    results = extract_best_poses(output_dir)
    
    print(f"\nTop 5 poses by affinity:")
    for i, r in enumerate(results[:5]):
        print(f"  {i+1}. {r['file']}: {r['affinity']:.3f} kcal/mol")
    
    # Save best pose as PDB
    if results:
        best = results[0]
        best_ligand_pdb = complex_dir / "best_ligand.pdb"
        with open(best_ligand_pdb, 'w') as f:
            f.write(pdbqt_to_pdb(best['content']))
        print(f"\nBest ligand pose saved to: {best_ligand_pdb}")
        
        # Create complex by combining protein and ligand
        protein_pdb = Path("pdb/protein_clean.pdb")
        complex_pdb = complex_dir / "protein_ligand_complex.pdb"
        
        with open(protein_pdb, 'r') as f:
            protein_content = f.read()
        
        # Remove END line from protein
        protein_lines = [l for l in protein_content.split('\n') if not l.startswith('END')]
        
        with open(complex_pdb, 'w') as f:
            f.write('\n'.join(protein_lines))
            f.write('\nTER\n')
            f.write(pdbqt_to_pdb(best['content']))
            f.write('\nEND\n')
        
        print(f"Complex PDB saved to: {complex_pdb}")
        
        # Summary
        print(f"\n=== Summary ===")
        print(f"Total docking results: {len(results)}")
        print(f"Best affinity: {best['affinity']:.3f} kcal/mol")
        print(f"Best pose from: {best['file']}")


if __name__ == "__main__":
    main()
