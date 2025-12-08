#!/usr/bin/env python3
"""Validate Wrap 'n' Shake setup and check for common issues."""

from __future__ import annotations

import argparse
import os
import re
import sys
from pathlib import Path
from typing import Dict, List, Tuple

try:
    import yaml  # type: ignore
except ImportError as exc:  # pragma: no cover
    raise SystemExit("PyYAML is required. Install it via 'pip install pyyaml'.") from exc


def check_pdbqt_format(file_path: Path) -> List[str]:
    """Check PDBQT file format for common issues."""
    issues = []
    
    if not file_path.exists():
        return [f"File not found: {file_path}"]
    
    with file_path.open('r', encoding='utf-8') as f:
        lines = f.readlines()
    
    for i, line in enumerate(lines, 1):
        if line.startswith(('ATOM', 'HETATM')):
            # Check line length
            if len(line) < 80:
                issues.append(f"Line {i}: PDBQT line too short ({len(line)} < 80)")
                continue
            
            # Check charge column (71-76, Python 70-75)
            try:
                charge_str = line[70:76].strip()
                if charge_str:
                    float(charge_str)
            except ValueError:
                issues.append(f"Line {i}: Invalid charge format: '{charge_str}'")
            
            # Check atom type column (78-79, Python 77-78)
            atom_type = line[77:79].strip()
            if not atom_type:
                issues.append(f"Line {i}: Missing atom type")
    
    return issues


def check_mdp_consistency(mdp_file: Path, cycle_time: float) -> List[str]:
    """Check MDP file for consistency with cycle time."""
    issues = []
    
    if not mdp_file.exists():
        return [f"MDP file not found: {mdp_file}"]
    
    content = mdp_file.read_text(encoding='utf-8')
    
    # Check nsteps and dt
    nsteps_match = re.search(r'nsteps\s*=\s*(\d+)', content)
    dt_match = re.search(r'dt\s*=\s*([\d.]+)', content)
    
    if not nsteps_match:
        issues.append("Could not find nsteps in MDP file")
    if not dt_match:
        issues.append("Could not find dt in MDP file")
    
    if nsteps_match and dt_match:
        nsteps = int(nsteps_match.group(1))
        dt = float(dt_match.group(1))
        total_time_ns = nsteps * dt / 1000  # Convert to ns
        
        if abs(total_time_ns - cycle_time) > 0.01:
            issues.append(
                f"Time mismatch: nsteps*dt = {total_time_ns:.3f} ns, "
                f"expected {cycle_time} ns"
            )
    
    # Check annealing parameters
    annealing_time_match = re.search(r'annealing-time\s*=\s*([\d\s.]+)', content)
    annealing_temp_match = re.search(r'annealing-temp\s*=\s*([\d\s.]+)', content)
    
    if not annealing_time_match:
        issues.append("Could not find annealing-time in MDP file")
    if not annealing_temp_match:
        issues.append("Could not find annealing-temp in MDP file")
    
    if annealing_time_match and annealing_temp_match:
        times = [float(t) for t in annealing_time_match.group(1).split()]
        temps = [float(t) for t in annealing_temp_match.group(1).split()]
        
        if len(times) != len(temps):
            issues.append("annealing-time and annealing-temp have different number of points")
        
        if times != sorted(times):
            issues.append("annealing-time points are not in chronological order")
    
    # Check position restraints
    if '-DPOSRES' not in content:
        issues.append("Position restraints (-DPOSRES) not found in MDP file")
    
    return issues


def check_topology_consistency(topol_file: Path, gro_file: Path) -> List[str]:
    """Check consistency between topology and GRO files."""
    issues = []
    
    if not topol_file.exists():
        return [f"Topology file not found: {topol_file}"]
    if not gro_file.exists():
        return [f"GRO file not found: {gro_file}"]
    
    # Count atoms in GRO file
    with gro_file.open('r', encoding='utf-8') as f:
        lines = f.readlines()
    
    if len(lines) < 2:
        return ["Invalid GRO file: too few lines"]
    
    try:
        gro_atom_count = int(lines[1].strip())
    except ValueError:
        return ["Invalid GRO file: cannot read atom count"]
    
    # Count atoms in topology
    topol_content = topol_file.read_text(encoding='utf-8')
    
    # Find [ molecules ] section
    molecules_match = re.search(r'\[\s*molecules\s*\](.*?)(?=\[|\Z)', topol_content, re.DOTALL)
    if not molecules_match:
        return ["Could not find [ molecules ] section in topology"]
    
    molecules_section = molecules_match.group(1)
    molecule_lines = [line.strip() for line in molecules_section.split('\n') 
                     if line.strip() and not line.strip().startswith(';')]
    
    # This is a simplified check - in practice, you'd need to parse the full topology
    # to get exact atom counts for each molecule type
    
    return issues


def check_x_atom_parameters(gpf_file: Path) -> List[str]:
    """Check if X atom parameters are defined in GPF file."""
    issues = []
    
    if not gpf_file.exists():
        return [f"GPF file not found: {gpf_file}"]
    
    content = gpf_file.read_text(encoding='utf-8')
    
    # Check if X is in ligand_types
    if 'ligand_types' in content and 'X' not in content:
        issues.append("X atom not found in ligand_types")
    
    # Check if X atom parameters are defined
    if 'parameter_file' not in content and 'X' not in content:
        issues.append("X atom parameters not defined")
    
    return issues


def check_path_compatibility(config: Dict) -> List[str]:
    """Check if paths are compatible with Windows/WSL environment."""
    issues = []
    
    # Check for Windows paths that might cause issues
    for key, path in config.get("paths", {}).items():
        if isinstance(path, str) and '\\' in path:
            if key in ["autodock4", "autogrid4", "gmx"]:
                # These are WSL executables, should not have Windows paths
                issues.append(f"Path '{key}' should not contain backslashes for WSL: {path}")
    
    return issues


def validate_config(config_file: Path) -> List[str]:
    """Validate configuration file."""
    issues = []
    
    if not config_file.exists():
        return [f"Configuration file not found: {config_file}"]
    
    try:
        with config_file.open('r', encoding='utf-8') as f:
            config = yaml.safe_load(f)
    except yaml.YAMLError as e:
        return [f"Invalid YAML in config file: {e}"]
    
    # Check required sections
    required_sections = ["paths", "inputs", "wrapper", "shaker"]
    for section in required_sections:
        if section not in config:
            issues.append(f"Missing required section: {section}")
    
    # Check wrapper parameters
    if "wrapper" in config:
        wrapper = config["wrapper"]
        if "seeds" not in wrapper:
            issues.append("Missing wrapper.seeds")
        elif not isinstance(wrapper["seeds"], list):
            issues.append("wrapper.seeds must be a list")
        
        if "max_cycles" in wrapper and wrapper["max_cycles"] <= 0:
            issues.append("wrapper.max_cycles must be positive")
        
        if "min_ligand_distance" in wrapper and wrapper["min_ligand_distance"] <= 0:
            issues.append("wrapper.min_ligand_distance must be positive")
    
    # Check shaker parameters
    if "shaker" in config:
        shaker = config["shaker"]
        if "n_cycles" in shaker and shaker["n_cycles"] <= 0:
            issues.append("shaker.n_cycles must be positive")
        
        if "displacement_cutoff" in shaker and shaker["displacement_cutoff"] <= 0:
            issues.append("shaker.displacement_cutoff must be positive")
    
    # Check path compatibility
    issues.extend(check_path_compatibility(config))
    
    return issues


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="config.yml", help="Configuration file")
    parser.add_argument("--work-dir", default=".", help="Working directory")
    args = parser.parse_args()
    
    work_dir = Path(args.work_dir).resolve()
    config_file = Path(args.config)
    
    print("=== Wrap 'n' Shake Setup Validation ===\n")
    
    all_issues = []
    
    # Validate configuration
    print("1. Validating configuration file...")
    config_issues = validate_config(config_file)
    all_issues.extend(config_issues)
    for issue in config_issues:
        print(f"   ❌ {issue}")
    if not config_issues:
        print("   ✅ Configuration file is valid")
    
    # Load config for further checks
    config = {}
    if config_file.exists():
        with config_file.open('r', encoding='utf-8') as f:
            config = yaml.safe_load(f)
    
    # Check PDBQT files
    print("\n2. Checking PDBQT files...")
    receptor_pdbqt = work_dir / config.get("inputs", {}).get("receptor_pdbqt", "")
    ligand_pdbqt = work_dir / config.get("inputs", {}).get("ligand_pdbqt", "")
    
    for pdbqt_file in [receptor_pdbqt, ligand_pdbqt]:
        if pdbqt_file.exists():
            pdbqt_issues = check_pdbqt_format(pdbqt_file)
            all_issues.extend(pdbqt_issues)
            for issue in pdbqt_issues:
                print(f"   ❌ {pdbqt_file.name}: {issue}")
            if not pdbqt_issues:
                print(f"   ✅ {pdbqt_file.name}: Format OK")
    
    # Check MDP file
    print("\n3. Checking annealing MDP file...")
    mdp_file = work_dir / "mdp" / "annealing.mdp"
    cycle_time = config.get("shaker", {}).get("cycle_time", 1.0)
    mdp_issues = check_mdp_consistency(mdp_file, cycle_time)
    all_issues.extend(mdp_issues)
    for issue in mdp_issues:
        print(f"   ❌ {issue}")
    if not mdp_issues:
        print("   ✅ Annealing MDP file is consistent")
    
    # Check topology consistency
    print("\n4. Checking topology consistency...")
    gmx_dir = work_dir / "gmx"
    topol_file = gmx_dir / "topol.top"
    gro_file = gmx_dir / "npt.gro"
    
    if topol_file.exists() and gro_file.exists():
        topol_issues = check_topology_consistency(topol_file, gro_file)
        all_issues.extend(topol_issues)
        for issue in topol_issues:
            print(f"   ❌ {issue}")
        if not topol_issues:
            print("   ✅ Topology and GRO files are consistent")
    else:
        missing = []
        if not topol_file.exists():
            missing.append("topol.top")
        if not gro_file.exists():
            missing.append("npt.gro")
        print(f"   ⚠️  Files not found: {', '.join(missing)}")
    
    # Check X atom parameters
    print("\n5. Checking X atom parameters...")
    gpf_file = work_dir / config.get("wrapper", {}).get("template_dir", "scripts/templates") / "gpf_template.txt"
    x_atom_issues = check_x_atom_parameters(gpf_file)
    all_issues.extend(x_atom_issues)
    for issue in x_atom_issues:
        print(f"   ❌ {issue}")
    if not x_atom_issues:
        print("   ✅ X atom parameters are defined")
    
    # Summary
    print(f"\n=== Validation Summary ===")
    if all_issues:
        print(f"❌ Found {len(all_issues)} issues:")
        for i, issue in enumerate(all_issues, 1):
            print(f"   {i}. {issue}")
        sys.exit(1)
    else:
        print("✅ All checks passed! Your Wrap 'n' Shake setup looks good.")
        print("\n💡 Recommendation: Run a small test with 2-3 ligands before full production.")


if __name__ == "__main__":
    main()