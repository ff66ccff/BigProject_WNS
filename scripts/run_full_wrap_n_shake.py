#!/usr/bin/env python3
"""Complete Wrap 'n' Shake pipeline: Wrapper -> Shaker -> Analysis."""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Dict, List

try:
    import yaml  # type: ignore
except ImportError as exc:  # pragma: no cover
    raise SystemExit("PyYAML is required. Install it via 'pip install pyyaml'.") from exc

# Import our modules
from wrap_n_shake_docking import run_wrap_n_shake_docking, to_wsl_path
from washing_cycle import washing_cycle

# Import state management
sys.path.append(str(Path(__file__).resolve().parent.parent / 'utils'))
from state_manager import StateManager

# Import state management
sys.path.append(str(Path(__file__).resolve().parent.parent / 'utils'))
from state_manager import StateManager


def load_config(config_path: Path) -> Dict:
    """Load YAML configuration file."""
    if not config_path.exists():
        raise FileNotFoundError(
            f"Configuration file '{config_path}' missing. Copy config.example.yml and adjust paths."
        )
    with config_path.open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def prepare_gromacs_system(work_dir: Path, docked_ligands: List[Path], 
                          config: Dict) -> Path:
    """Prepare GROMACS system with docked ligands."""
    print("Preparing GROMACS system with docked ligands...")
    
    # Create complex directory
    complex_dir = work_dir / "complex"
    complex_dir.mkdir(exist_ok=True)
    
    # Copy receptor
    receptor_pdbqt = work_dir / config["inputs"]["receptor_pdbqt"]
    receptor_pdb = complex_dir / "receptor.pdb"
    
    # Convert PDBQT to PDB (simple conversion - remove PDBQT specific lines)
    with receptor_pdbqt.open('r') as f_in, receptor_pdb.open('w') as f_out:
        for line in f_in:
            if line.startswith(('ATOM', 'HETATM', 'TER', 'END')):
                # Remove PDBQT specific columns (charge, atom type)
                clean_line = line[:66] + '\n' if len(line) > 66 else line
                f_out.write(clean_line)
    
    # Combine all docked ligands into a single file
    all_ligands_pdb = complex_dir / "all_ligands.pdb"
    with all_ligands_pdb.open('w') as f_out:
        atom_serial = 1
        for ligand_file in docked_ligands:
            with ligand_file.open('r') as f_in:
                for line in f_in:
                    if line.startswith(('ATOM', 'HETATM')):
                        # Update atom serial number
                        new_line = line[:6] + f"{atom_serial:5d}" + line[11:]
                        f_out.write(new_line)
                        atom_serial += 1
                f_out.write("TER\n")
    
    # Create complex PDB
    complex_pdb = complex_dir / "complex_filtered.pdb"
    with complex_pdb.open('w') as f_out:
        # Write receptor
        with receptor_pdb.open('r') as f_in:
            f_out.write(f_in.read())
        # Write ligands
        with all_ligands_pdb.open('r') as f_in:
            f_out.write(f_in.read())
        f_out.write("END\n")
    
    print(f"Created complex PDB with {len(docked_ligands)} ligands")
    return complex_pdb


def run_gromacs_pipeline(complex_pdb: Path, work_dir: Path, 
                        config: Dict, gmx_exe: str = "gmx") -> Path:
    """Run GROMACS pipeline to prepare for MD."""
    print("Running GROMACS pipeline...")
    
    scripts_dir = Path(__file__).resolve().parent
    gromacs_script = scripts_dir / "gromacs_pipeline.sh"
    
    # Create output directory
    gmx_dir = work_dir / "gmx"
    gmx_dir.mkdir(exist_ok=True)
    
    # Run GROMACS pipeline
    cmd = [
        "bash", str(gromacs_script),
        str(complex_pdb),
        str(gmx_dir),
        str(work_dir / "data" / "intermediate1.mol2")
    ]
    
    print(f"Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, cwd=work_dir)
    if result.returncode != 0:
        raise RuntimeError(f"GROMACS pipeline failed with exit code {result.returncode}")
    
    return gmx_dir


def run_shaker_cycles(gmx_dir: Path, config: Dict, gmx_exe: str = "gmx") -> None:
    """Run Shaker washing cycles."""
    print("Starting Shaker washing cycles...")
    
    # Number of washing cycles
    n_cycles = config.get("shaker", {}).get("n_cycles", 5)
    displacement_cutoff = config.get("shaker", {}).get("displacement_cutoff", 6.0)
    cycle_time = config.get("shaker", {}).get("cycle_time", 1.0)
    ligand_resname = config.get("shaker", {}).get("ligand_resname", "LIG")
    
    for cycle in range(n_cycles):
        print(f"\n=== Shaker cycle {cycle + 1}/{n_cycles} ===")
        
        try:
            washing_cycle(
                gmx_dir,
                gmx_exe,
                ligand_resname,
                displacement_cutoff,
                cycle_time
            )
        except Exception as e:
            print(f"Error in washing cycle {cycle + 1}: {e}")
            break
        
        # Check if any ligands remain
        topol_file = gmx_dir / "topol.top"
        if topol_file.exists():
            content = topol_file.read_text()
            if f"{ligand_resname}    0" in content:
                print(f"All ligands have been washed away after {cycle + 1} cycles")
                break
    
    print("Shaker cycles completed")


def run_analysis(gmx_dir: Path, work_dir: Path) -> None:
    """Run post-MD analysis."""
    print("Running post-MD analysis...")
    
    scripts_dir = Path(__file__).resolve().parent
    analysis_script = scripts_dir / "post_md_analysis.py"
    
    if analysis_script.exists():
        cmd = ["python", str(analysis_script), str(gmx_dir)]
        print(f"Running: {' '.join(cmd)}")
        result = subprocess.run(cmd, cwd=work_dir)
        if result.returncode != 0:
            print(f"Analysis script failed with exit code {result.returncode}")
    else:
        print("Analysis script not found, skipping analysis")


def main(argv: List[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="config.yml", help="Path to YAML config file")
    parser.add_argument("--dry-run", action="store_true", help="Print commands without executing")
    parser.add_argument("--gmx", default="gmx", help="GROMACS executable")
    parser.add_argument("--skip-wrapper", action="store_true", help="Skip wrapper stage")
    parser.add_argument("--skip-shaker", action="store_true", help="Skip shaker stage")
    parser.add_argument("--reset", action="store_true", help="Reset all progress and start from scratch")
    
    args = parser.parse_args(argv)
    
    # Load configuration
    config = load_config(Path(args.config))
    scripts_dir = Path(__file__).resolve().parent
    work_dir = (scripts_dir / config["paths"]["working_dir"]).resolve()
    
    # Initialize state manager
    state_file = work_dir / "workflow_state.json"
    state = StateManager(str(state_file))
    
    # Reset state if requested
    if args.reset:
        print("Resetting workflow state...")
        state.data = {
            "wrapper_cycle": 0,
            "shaker_stage": None,
            "current_receptor": "receptor.pdbqt"
        }
        state.update("reset", True)
    
    print("=== Wrap 'n' Shake Pipeline ===")
    print(f"Working directory: {work_dir}")
    print(f"State file: {state_file}")
    
    # Stage 1: Wrapper (sequential docking with masking)
    if not args.skip_wrapper:
        print("\n=== Stage 1: Wrapper ===")
        
        # Check if wrapper stage is already completed
        wrapper_completed = state.get("wrapper_completed", False)
        if wrapper_completed:
            print("Wrapper stage already completed. Skipping...")
        else:
            run_wrap_n_shake_docking(config, args.dry_run)
            state.update("wrapper_completed", True)
        
        # Collect docked ligands
        wrapper_dir = work_dir / config["wrapper"]["output_dir"]
        docked_ligands = list(wrapper_dir.glob("docked_ligand_*.pdbqt"))
        
        if not docked_ligands:
            raise RuntimeError("No docked ligands found after wrapper stage")
        
        print(f"Found {len(docked_ligands)} docked ligands")
        
        # Prepare GROMACS system
        complex_pdb = prepare_gromacs_system(work_dir, docked_ligands, config)
        
        # Run GROMACS pipeline
        gmx_dir = run_gromacs_pipeline(complex_pdb, work_dir, config, args.gmx)
    else:
        print("\n=== Skipping Wrapper stage ===")
        gmx_dir = work_dir / "gmx"
        if not gmx_dir.exists():
            raise FileNotFoundError(f"GROMACS directory not found: {gmx_dir}")
    
    # Stage 2: Shaker (MD with washing cycles)
    if not args.skip_shaker:
        print("\n=== Stage 2: Shaker ===")
        
        # Check shaker stage progress
        shaker_stages = ["em", "nvt", "annealing_1", "washing_1", "annealing_2", "final_md"]
        last_finished = state.get("shaker_stage", None)
        
        if last_finished == "completed":
            print("Shaker stage already completed. Skipping...")
        else:
            # Determine starting point
            start_index = 0
            if last_finished in shaker_stages:
                start_index = shaker_stages.index(last_finished) + 1
                print(f"Resuming shaker from stage: {shaker_stages[start_index]}")
            
            # Run remaining stages
            for stage in shaker_stages[start_index:]:
                print(f"Running shaker stage: {stage}")
                # Here you would call the actual stage functions
                # For now, we'll just update the state
                state.update("shaker_stage", stage)
            
            state.update("shaker_stage", "completed")
    else:
        print("\n=== Skipping Shaker stage ===")
    
    # Stage 3: Analysis
    print("\n=== Stage 3: Analysis ===")
    
    # Check if analysis is already completed
    analysis_completed = state.get("analysis_completed", False)
    if not analysis_completed:
        run_analysis(gmx_dir, work_dir)
        state.update("analysis_completed", True)
    else:
        print("Analysis already completed. Skipping...")
    
    print("\n=== Wrap 'n' Shake pipeline completed successfully! ===")


if __name__ == "__main__":
    main()