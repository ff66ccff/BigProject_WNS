#!/usr/bin/env python3
"""Wrap 'n' Shake docking pipeline: sequential docking with atom masking."""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path
from string import Template
from typing import Dict, List

try:
    import yaml  # type: ignore
except ImportError as exc:  # pragma: no cover
    raise SystemExit("PyYAML is required. Install it via 'pip install pyyaml'.") from exc

# Import our masking function
from mask_pdbqt import mask_receptor


def to_wsl_path(path: Path) -> str:
    """Convert Windows path to WSL path."""
    drive = path.drive.rstrip(":")
    if not drive:
        return str(path.as_posix())
    return "/mnt/" + drive.lower() + path.as_posix()[2:]


def load_config(config_path: Path) -> Dict:
    """Load YAML configuration file."""
    if not config_path.exists():
        raise FileNotFoundError(
            f"Configuration file '{config_path}' missing. Copy config.example.yml and adjust paths."
        )
    with config_path.open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def render_template(template_path: Path, mapping: Dict[str, str]) -> str:
    """Render template file with given mapping."""
    template = Template(template_path.read_text(encoding="utf-8"))
    return template.safe_substitute(mapping)


def write_file(path: Path, content: str) -> None:
    """Write content to file, creating parent directories if needed."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")


def windows_command_exists(cmd_name: str) -> bool:
    """Check whether a command exists on Windows PATH."""
    import shutil
    return shutil.which(cmd_name) is not None


def wsl_command_exists(cmd_name: str) -> bool:
    """Check whether a command exists in WSL PATH."""
    try:
        result = subprocess.run(
            ["wsl", "which", cmd_name],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        return result.returncode == 0
    except Exception:
        return False


def build_command(executable: str, args: List[str], use_wsl: bool = False) -> List[str]:
    """Build command list, wrapping with wsl if needed."""
    if executable.startswith("/") or use_wsl:
        return ["wsl", executable, *args]
    return [executable, *args]


def run_command(cmd: List[str], dry_run: bool, cwd: Path | None = None) -> None:
    """Run command and handle errors."""
    print(" ".join(cmd))
    if dry_run:
        return
    result = subprocess.run(cmd, check=False, cwd=cwd)
    if result.returncode != 0:
        raise RuntimeError(f"Command '{cmd}' failed with exit code {result.returncode}")


def extract_best_pose(dlg_file: Path, output_pdbqt: Path) -> None:
    """Extract the best (lowest energy) pose from AutoDock DLG file."""
    best_pose_lines = []
    found_best = False
    
    with dlg_file.open('r', encoding='utf-8') as f:
        lines = f.readlines()
    
    # Find the best pose section
    for i, line in enumerate(lines):
        if "DOCKED: USER    Run = 0" in line:
            # Found the best pose, extract coordinates
            j = i
            while j < len(lines) and not lines[j].startswith("DOCKED: USER    Run = 1"):
                if lines[j].startswith("DOCKED: ATOM"):
                    best_pose_lines.append(lines[j][8:])  # Remove "DOCKED: " prefix
                j += 1
            found_best = True
            break
    
    if not found_best:
        raise RuntimeError(f"Could not find best pose in {dlg_file}")
    
    # Write the best pose to PDBQT file
    with output_pdbqt.open('w', encoding='utf-8') as f:
        # Write PDBQT header
        f.write("REMARK  Best pose from AutoDock\n")
        f.write("TORSDOF 0\n")
        # Write atoms
        for line in best_pose_lines:
            f.write(line)
    
    print(f"Extracted best pose to {output_pdbqt}")


def run_wrap_n_shake_docking(config: Dict, dry_run: bool = False) -> None:
    """Run Wrap 'n' Shake docking pipeline."""
    scripts_dir = Path(__file__).resolve().parent
    working_dir = (scripts_dir / config["paths"]["working_dir"]).resolve()
    template_dir = (working_dir / config["wrapper"]["template_dir"]).resolve()
    
    receptor_pdbqt = working_dir / config["inputs"]["receptor_pdbqt"]
    ligand_pdbqt = working_dir / config["inputs"]["ligand_pdbqt"]
    output_dir = working_dir / config["wrapper"]["output_dir"]
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Copy original receptor to working directory
    receptor_current = output_dir / "receptor_current.pdbqt"
    shutil.copy2(receptor_pdbqt, receptor_current)
    
    # Setup AutoGrid
    gpf_path = output_dir / "autogrid.gpf"
    gridfld = output_dir / "protein.maps.fld"
    
    # Generate map definitions for each ligand type
    ligand_types = config["inputs"]["ligand_types"].split()
    map_definitions = []
    for ligand_type in ligand_types:
        map_definitions.append(f"map {ligand_type}.map")
    
    gpf_mapping = {
        "npts_x": str(config["autogrid"]["npts"][0]),
        "npts_y": str(config["autogrid"]["npts"][1]),
        "npts_z": str(config["autogrid"]["npts"][2]),
        "gridfld": gridfld.name,
        "receptor": receptor_current.name,
        "center_x": str(config["autogrid"]["center"][0]),
        "center_y": str(config["autogrid"]["center"][1]),
        "center_z": str(config["autogrid"]["center"][2]),
        "spacing": str(config["autogrid"]["spacing"]),
        "ligand_types": config["inputs"]["ligand_types"],
        "map_definitions": "\n".join(map_definitions),
    }
    
    gpf_content = render_template(template_dir / "gpf_template.txt", gpf_mapping)
    write_file(gpf_path, gpf_content)
    
    # Run AutoGrid
    autogrid_exe = config["paths"]["autogrid4"]
    use_wsl_autogrid = not windows_command_exists(autogrid_exe) and wsl_command_exists(autogrid_exe)
    
    if use_wsl_autogrid:
        wsl_output_dir = to_wsl_path(output_dir)
        autogrid_cmd = ["wsl", "bash", "-c", f"cd {wsl_output_dir} && {autogrid_exe} -p {gpf_path.name} -l autogrid.log"]
    else:
        autogrid_cmd = build_command(autogrid_exe, ["-p", gpf_path.name, "-l", "autogrid.log"], use_wsl=False)
    
    run_command(autogrid_cmd, dry_run, cwd=output_dir if not use_wsl_autogrid else None)
    
    # Sequential docking loop
    seeds = config["wrapper"]["seeds"]
    docked_ligands = []  # Store paths to docked ligand files
    
    for i, seed in enumerate(seeds):
        print(f"\n=== Docking cycle {i+1}/{len(seeds)} (seed: {seed}) ===")
        
        # Generate DPF file
        dpf_template = template_dir / "dpf_template.txt"
        map_definitions = "\n".join([
            f"map {ligand_type}.map" for ligand_type in config["inputs"]["ligand_types"].split()
        ])
        
        mapping = {
            "seed": str(seed),
            "ligand_types": config["inputs"]["ligand_types"],
            "gridfld": gridfld.name,
            "maps": map_definitions,
            "ligand_pdbqt": ligand_pdbqt.name,
            "receptor": receptor_current.stem,
            "center_x": str(config["autogrid"]["center"][0]),
            "center_y": str(config["autogrid"]["center"][1]),
            "center_z": str(config["autogrid"]["center"][2]),
        }
        
        dpf_content = render_template(dpf_template, mapping)
        dpf_path = output_dir / f"wrapper_{seed}.dpf"
        write_file(dpf_path, dpf_content)
        
        # Run AutoDock
        autodock_exe = config["paths"]["autodock4"]
        use_wsl_autodock = not windows_command_exists(autodock_exe) and wsl_command_exists(autodock_exe)
        
        dlg_path = output_dir / f"wrapper_{seed}.dlg"
        docked_ligand_path = output_dir / f"docked_ligand_{seed}.pdbqt"
        
        if use_wsl_autodock:
            wsl_output_dir = to_wsl_path(output_dir)
            cmd = ["wsl", "bash", "-c", f"cd {wsl_output_dir} && {autodock_exe} -p {dpf_path.name} -l {dlg_path.name}"]
        else:
            cmd = build_command(autodock_exe, ["-p", dpf_path.name, "-l", dlg_path.name], use_wsl=False)
        
        run_command(cmd, dry_run, cwd=output_dir if not use_wsl_autodock else None)
        
        # Extract best pose
        extract_best_pose(dlg_path, docked_ligand_path)
        docked_ligands.append(docked_ligand_path)
        
        # Mask receptor atoms for next iteration
        if i < len(seeds) - 1:  # Don't mask after the last docking
            print(f"Masking receptor atoms near docked ligand {seed}...")
            receptor_next = output_dir / f"receptor_masked_{seed}.pdbqt"
            mask_receptor(receptor_current, [docked_ligand_path], receptor_next, cutoff=3.5)
            receptor_current = receptor_next  # Use masked receptor for next iteration
    
    print(f"\n=== Wrap 'n' Shake docking completed ===")
    print(f"Docked {len(docked_ligands)} ligands:")
    for ligand_path in docked_ligands:
        print(f"  - {ligand_path}")


def main(argv: List[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="config.yml", help="Path to YAML config file")
    parser.add_argument("--dry-run", action="store_true", help="Print commands without executing")
    args = parser.parse_args(argv)
    
    config = load_config(Path(args.config))
    run_wrap_n_shake_docking(config, args.dry_run)


if __name__ == "__main__":
    main()