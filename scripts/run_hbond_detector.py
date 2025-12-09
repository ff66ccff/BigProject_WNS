#!/usr/bin/env python3
"""Hydrogen Bond Detector: Complete WnS pipeline for hydrogen bond detection."""

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
from analyze_hbonds import calculate_wns_score, save_results_to_csv
from washing_cycle import washing_cycle


def load_config(config_path: Path) -> Dict:
    """Load YAML configuration file."""
    if not config_path.exists():
        raise FileNotFoundError(
            f"Configuration file '{config_path}' missing. Copy config.example.yml and adjust paths."
        )
    with config_path.open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def run_command(cmd: List[str], cwd: Path | None = None, dry_run: bool = False) -> None:
    """Run command and handle errors."""
    print(f"Running: {' '.join(cmd)}")
    if dry_run:
        return
    
    result = subprocess.run(cmd, cwd=cwd, check=False)
    if result.returncode != 0:
        raise RuntimeError(f"Command '{cmd}' failed with exit code {result.returncode}")


def run_wrapper_stage(config: Dict, dry_run: bool = False) -> Path:
    """Run Wrapper stage to generate ligand poses on protein surface."""
    print("\n=== Stage 1: Wrapper - Generating ligand poses ===")
    
    scripts_dir = Path(__file__).resolve().parent
    working_dir = (scripts_dir / config["paths"]["working_dir"]).resolve()
    
    # Run AutoDock batch
    wrapper_script = scripts_dir / "run_autodock_batch.py"
    cmd = ["python", str(wrapper_script), "--config", str(scripts_dir / "config.yml")]
    if dry_run:
        cmd.append("--dry-run")
    
    run_command(cmd, cwd=scripts_dir, dry_run=dry_run)
    
    # Check for output
    output_dir = working_dir / config["wrapper"]["output_dir"]
    wrapped_complex = output_dir / "wrapped_complex.pdb"
    
    if not dry_run and not wrapped_complex.exists():
        raise RuntimeError(f"Wrapper stage failed: {wrapped_complex} not found")
    
    print(f"Wrapper stage completed: {wrapped_complex}")
    return wrapped_complex


def run_shaker_stage(config: Dict, wrapped_complex: Path, dry_run: bool = False) -> Path:
    """Run Shaker stage with MD simulation to test hydrogen bond stability."""
    print("\n=== Stage 2: Shaker - Testing hydrogen bond stability ===")
    
    scripts_dir = Path(__file__).resolve().parent
    working_dir = (scripts_dir / config["paths"]["working_dir"]).resolve()
    gmx_dir = working_dir / config["gromacs"]["workdir"]
    
    # Copy wrapped complex to GROMACS directory
    gmx_complex = gmx_dir / "complex.pdb"
    if not dry_run:
        gmx_dir.mkdir(parents=True, exist_ok=True)
        shutil.copy2(wrapped_complex, gmx_complex)
    
    # Run GROMACS pipeline (simplified version)
    gmx_script = scripts_dir / "gromacs_pipeline.py"
    if gmx_script.exists():
        cmd = ["python", str(gmx_script), "--input", str(gmx_complex), "--config", str(scripts_dir / "config.yml")]
        if dry_run:
            cmd.append("--dry-run")
        run_command(cmd, cwd=scripts_dir, dry_run=dry_run)
    else:
        print("Warning: gromacs_pipeline.py not found, skipping GROMACS setup")
    
    # Run washing cycle
    if not dry_run:
        washing_cycle(
            gmx_dir,
            gmx_exe=config["paths"]["gmx"],
            ligand_resname="LIG",
            displacement_cutoff=6.0,
            cycle_time=1.0
        )
    
    # Check for output
    final_complex = gmx_dir / "npt.gro"
    if not dry_run and not final_complex.exists():
        raise RuntimeError(f"Shaker stage failed: {final_complex} not found")
    
    # Convert GRO to PDB for analysis
    final_pdb = gmx_dir / "final_complex.pdb"
    if not dry_run:
        # Use editconf to convert GRO to PDB
        cmd = [config["paths"]["gmx"], "editconf", "-f", "npt.gro", "-o", "final_complex.pdb"]
        run_command(cmd, cwd=gmx_dir)
    
    print(f"Shaker stage completed: {final_pdb}")
    return final_pdb


def run_analysis_stage(config: Dict, final_complex: Path, dry_run: bool = False) -> None:
    """Run Analysis stage to count hydrogen bonds and generate visualization."""
    print("\n=== Stage 3: Analysis - Counting hydrogen bonds and generating visualization ===")
    
    scripts_dir = Path(__file__).resolve().parent
    working_dir = (scripts_dir / config["paths"]["working_dir"]).resolve()
    
    # Run hydrogen bond analysis
    analysis_script = scripts_dir / "analyze_hbonds.py"
    output_prefix = working_dir / "hbond_results"
    
    cmd = [
        "python", str(analysis_script), str(final_complex),
        "-o", str(output_prefix),
        "--csv"
    ]
    
    if dry_run:
        cmd.append("--dry-run")
    
    run_command(cmd, cwd=working_dir, dry_run=dry_run)
    
    # Check for outputs
    csv_file = working_dir / "hbond_results.csv"
    
    if not dry_run:
        if csv_file.exists():
            print(f"Analysis results: {csv_file}")
            print(f"Only ligands with hydrogen bonds are included in the output")
    
    print("Analysis stage completed")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="config.yml", help="Path to YAML config file")
    parser.add_argument("--dry-run", action="store_true", help="Print commands without executing")
    parser.add_argument("--stage", choices=["wrapper", "shaker", "analysis", "all"], 
                       default="all", help="Which stage to run")
    parser.add_argument("--input", type=Path, help="Input PDB file (skip wrapper stage)")
    
    args = parser.parse_args()
    
    # Load configuration
    scripts_dir = Path(__file__).resolve().parent
    config_path = scripts_dir / args.config
    config = load_config(config_path)
    
    try:
        if args.stage in ["wrapper", "all"]:
            wrapped_complex = run_wrapper_stage(config, args.dry_run)
        else:
            if args.input:
                wrapped_complex = args.input
            else:
                # Use default output from wrapper stage
                working_dir = (scripts_dir / config["paths"]["working_dir"]).resolve()
                output_dir = working_dir / config["wrapper"]["output_dir"]
                wrapped_complex = output_dir / "wrapped_complex.pdb"
                
                if not wrapped_complex.exists():
                    raise FileNotFoundError(f"Input file not found: {wrapped_complex}")
        
        if args.stage in ["shaker", "all"]:
            final_complex = run_shaker_stage(config, wrapped_complex, args.dry_run)
        else:
            if args.input:
                final_complex = args.input
            else:
                # Use default output from shaker stage
                working_dir = (scripts_dir / config["paths"]["working_dir"]).resolve()
                gmx_dir = working_dir / config["gromacs"]["workdir"]
                final_complex = gmx_dir / "final_complex.pdb"
                
                if not final_complex.exists():
                    raise FileNotFoundError(f"Input file not found: {final_complex}")
        
        if args.stage in ["analysis", "all"]:
            run_analysis_stage(config, final_complex, args.dry_run)
        
        print("\n=== Hydrogen Bond Detector pipeline completed successfully ===")
        
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()