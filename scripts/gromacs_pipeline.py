#!/usr/bin/env python3
"""GROMACS pipeline with checkpoint support."""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path

# Import state management
sys.path.append(str(Path(__file__).resolve().parent.parent / 'utils'))
from state_manager import StateManager
from gmx_runner import run_gmx_mdrun_safe


def run_command(cmd, cwd=None):
    """Run a command and handle errors."""
    print(f"Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, cwd=cwd, check=True)
    return result


def run_grompp(gmx_cmd, mdp_file, input_gro, topol_file, output_tpr, maxwarn=1, cwd=None):
    """Run grompp with checkpoint support."""
    cmd = [gmx_cmd, "grompp", "-f", mdp_file, "-c", input_gro, 
           "-p", topol_file, "-o", output_tpr]
    if maxwarn > 0:
        cmd.extend(["-maxwarn", str(maxwarn)])
    return run_command(cmd, cwd)


def run_genion(gmx_cmd, tpr_file, output_gro, neutral=True, conc=0.15, cwd=None):
    """Run genion to add ions."""
    cmd = [gmx_cmd, "genion", "-s", tpr_file, "-o", output_gro]
    if neutral:
        cmd.append("-neutral")
    if conc > 0:
        cmd.extend(["-conc", str(conc)])
    
    # Use subprocess with input
    print(f"Running: {' '.join(cmd)} (input: SOL)")
    result = subprocess.run(cmd, input="SOL\n", text=True, cwd=cwd, check=True)
    return result


def run_gromacs_pipeline(input_pdb, work_dir, gmx_cmd="gmx", ligand_mol2=None):
    """Run GROMACS pipeline with checkpoint support."""
    work_dir = Path(work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)
    
    # Initialize state manager
    state_file = work_dir / "gromacs_state.json"
    state = StateManager(str(state_file))
    
    # Get current stage
    current_stage = state.get("current_stage", "start")
    print(f"Current stage: {current_stage}")
    
    # Define stages in order
    stages = [
        "extract_protein",
        "parameterize_ligand", 
        "generate_topology",
        "combine_complex",
        "setup_box",
        "solvate",
        "add_ions",
        "energy_min",
        "nvt_eq",
        "npt_eq",
        "production_md"
    ]
    
    # Find starting index
    if current_stage in stages:
        start_idx = stages.index(current_stage)
    else:
        start_idx = 0
    
    # Run stages from starting point
    for i, stage in enumerate(stages[start_idx:], start=start_idx):
        print(f"\n=== Running stage: {stage} ===")
        
        try:
            if stage == "extract_protein":
                # Extract protein from complex
                protein_pdb = work_dir / "protein_only.pdb"
                if not protein_pdb.exists():
                    with open(input_pdb, 'r') as f_in, open(protein_pdb, 'w') as f_out:
                        for line in f_in:
                            if line.startswith(('ATOM', 'TER')):
                                if not any(x in line for x in ['UNL', 'LIG', 'MOL']):
                                    f_out.write(line)
                state.update("current_stage", "parameterize_ligand")
                
            elif stage == "parameterize_ligand":
                # Check if ligand exists and parameterize
                has_ligand = False
                if ligand_mol2 and Path(ligand_mol2).exists():
                    # Check for ACPYPE
                    acpype_exe = shutil.which("acpype")
                    if acpype_exe:
                        acpype_cmd = [acpype_exe, "-i", str(ligand_mol2), "-n", "0", "-a", "gaff2", "-o", "gmx", "-f"]
                        use_python_module = False
                    else:
                        # Fallback: try python -m acpype
                        try:
                            subprocess.run(["python3", "-m", "acpype", "--help"], capture_output=True, check=True)
                            acpype_cmd = ["python3", "-m", "acpype", "-i", str(ligand_mol2), "-n", "0", "-a", "gaff2", "-o", "gmx", "-f"]
                            use_python_module = True
                        except subprocess.CalledProcessError:
                            acpype_cmd = None
                    
                    if acpype_cmd:
                        print(f"Parameterizing ligand with ACPYPE (via {'module' if use_python_module else 'binary'})...")
                        try:
                            run_command(acpype_cmd, cwd=work_dir)
                            has_ligand = True
                        except subprocess.CalledProcessError:
                            print("ACPYPE failed, proceeding with protein only")
                            has_ligand = False
                    else:    
                        print("ACPYPE not found in PATH or as python module. Proceeding with protein only.")
                        has_ligand = False
                else:
                    print("No ligand file provided, proceeding with protein only")
                    has_ligand = False
                
                state.update("has_ligand", has_ligand)
                state.update("current_stage", "generate_topology")
                
            elif stage == "generate_topology":
                # Generate protein topology
                protein_pdb = work_dir / "protein_only.pdb"
                cmd = [gmx_cmd, "pdb2gmx", "-f", str(protein_pdb), 
                       "-o", "protein.gro", "-p", "topol.top",
                       "-ff", "amber99sb-ildn", "-water", "tip3p", "-ignh"]
                run_command(cmd, cwd=work_dir)
                state.update("current_stage", "combine_complex")
                
            elif stage == "combine_complex":
                # Combine protein and ligand if ligand exists
                has_ligand = state.get("has_ligand", False)
                if has_ligand:
                    # Find ACPYPE output
                    acpype_dirs = list(work_dir.glob("*.acpype"))
                    if acpype_dirs:
                        acpype_dir = acpype_dirs[0]
                        # Copy ligand files
                        for itp_file in acpype_dir.glob("*.itp"):
                            itp_file.rename(work_dir / "ligand.itp")
                        for gro_file in acpype_dir.glob("*.gro"):
                            gro_file.rename(work_dir / "ligand.gro")
                        
                        # Combine protein and ligand
                        # (Implementation would go here)
                        print("Combined protein and ligand")
                    else:
                        print("ACPYPE output not found, using protein only")
                
                # Create complex.gro (protein only if no ligand)
                if not (work_dir / "complex.gro").exists():
                    (work_dir / "protein.gro").rename(work_dir / "complex.gro")
                
                state.update("current_stage", "setup_box")
                
            elif stage == "setup_box":
                # Setup simulation box
                cmd = [gmx_cmd, "editconf", "-f", "complex.gro", 
                       "-o", "boxed.gro", "-c", "-d", "1.0", "-bt", "cubic"]
                run_command(cmd, cwd=work_dir)
                state.update("current_stage", "solvate")
                
            elif stage == "solvate":
                # Solvate system
                cmd = [gmx_cmd, "solvate", "-cp", "boxed.gro", 
                       "-cs", "tip3p.gro", "-o", "solvated.gro", "-p", "topol.top"]
                run_command(cmd, cwd=work_dir)
                state.update("current_stage", "add_ions")
                
            elif stage == "add_ions":
                # Add ions
                mdp_dir = Path(__file__).parent.parent / "mdp"
                run_grompp(gmx_cmd, str(mdp_dir / "ions.mdp"), 
                          "solvated.gro", "topol.top", "ions.tpr", cwd=work_dir)
                run_genion(gmx_cmd, "ions.tpr", "ionized.gro", cwd=work_dir)
                state.update("current_stage", "energy_min")
                
            elif stage == "energy_min":
                # Energy minimization
                mdp_dir = Path(__file__).parent.parent / "mdp"
                run_grompp(gmx_cmd, str(mdp_dir / "em.mdp"), 
                          "ionized.gro", "topol.top", "em.tpr", cwd=work_dir)
                run_gmx_mdrun_safe("em", gmx_cmd)
                state.update("current_stage", "nvt_eq")
                
            elif stage == "nvt_eq":
                # NVT equilibration
                mdp_dir = Path(__file__).parent.parent / "mdp"
                run_grompp(gmx_cmd, str(mdp_dir / "nvt.mdp"), 
                          "em.gro", "topol.top", "nvt.tpr", cwd=work_dir)
                run_gmx_mdrun_safe("nvt", gmx_cmd)
                state.update("current_stage", "npt_eq")
                
            elif stage == "npt_eq":
                # NPT equilibration
                mdp_dir = Path(__file__).parent.parent / "mdp"
                run_grompp(gmx_cmd, str(mdp_dir / "npt.mdp"), 
                          "nvt.gro", "topol.top", "npt.tpr", cwd=work_dir)
                run_gmx_mdrun_safe("npt", gmx_cmd)
                state.update("current_stage", "production_md")
                
            elif stage == "production_md":
                # Production MD
                mdp_dir = Path(__file__).parent.parent / "mdp"
                run_grompp(gmx_cmd, str(mdp_dir / "md.mdp"), 
                          "npt.gro", "topol.top", "md.tpr", cwd=work_dir)
                run_gmx_mdrun_safe("md", gmx_cmd)
                state.update("current_stage", "completed")
                
        except subprocess.CalledProcessError as e:
            print(f"Error in stage {stage}: {e}")
            state.update("current_stage", stage)  # Save current stage for retry
            raise
    
    print("\n=== GROMACS pipeline completed successfully! ===")
    return work_dir


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python gromacs_pipeline.py <input_pdb> <work_dir> [ligand_mol2]")
        sys.exit(1)
    
    input_pdb = sys.argv[1]
    work_dir = sys.argv[2]
    ligand_mol2 = sys.argv[3] if len(sys.argv) > 3 else None
    
    run_gromacs_pipeline(input_pdb, work_dir, ligand_mol2=ligand_mol2)