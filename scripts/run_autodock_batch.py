#!/usr/bin/env python3
"""Generate AutoGrid/AutoDock input files and run batch docking jobs."""

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

REPO_ROOT = Path(__file__).resolve().parents[1]
CONFIG_PATH = REPO_ROOT / "scripts" / "config.yml"


def to_wsl_path(path: Path) -> str:
    drive = path.drive.rstrip(":")
    if not drive:
        return str(path.as_posix())
    return "/mnt/" + drive.lower() + path.as_posix()[2:]


def load_config(config_path: Path) -> Dict:
    if not config_path.exists():
        raise FileNotFoundError(
            f"Configuration file '{config_path}' missing. Copy config.example.yml and adjust paths."
        )
    with config_path.open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def render_template(template_path: Path, mapping: Dict[str, str]) -> str:
    template = Template(template_path.read_text(encoding="utf-8"))
    return template.safe_substitute(mapping)


def write_file(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")


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


def windows_command_exists(cmd_name: str) -> bool:
    """Check whether a command exists on Windows PATH."""
    import shutil
    return shutil.which(cmd_name) is not None


def build_command(executable: str, args: List[str], use_wsl: bool = False) -> List[str]:
    """Build command list, wrapping with wsl if needed.
    
    Args:
        executable: The executable name or path.
        args: Command arguments.
        use_wsl: If True, force WSL execution.
    """
    if executable.startswith("/") or use_wsl:
        return ["wsl", executable, *args]
    return [executable, *args]


def run_command(cmd: List[str], dry_run: bool, cwd: Path | None = None) -> None:
    print(" ".join(cmd))
    if dry_run:
        return
    result = subprocess.run(cmd, check=False, cwd=cwd)
    if result.returncode != 0:
        raise RuntimeError(f"Command '{cmd}' failed with exit code {result.returncode}")


def merge_complex_files(receptor_pdbqt: Path, ligand_pdbqt_files: List[Path], output_pdb: Path) -> None:
    """Merge receptor and ligand PDBQT files into a single PDB file."""
    print(f"Merging complex files to {output_pdb}")
    
    with output_pdb.open('w', encoding='utf-8') as out_file:
        # Copy receptor atoms (convert from PDBQT to PDB format)
        with receptor_pdbqt.open('r', encoding='utf-8') as receptor_file:
            for line in receptor_file:
                if line.startswith(('ATOM', 'HETATM')):
                    # Convert PDBQT to PDB format (remove charge and PDBQT specific fields)
                    pdb_line = line[:66] + '\n'  # Truncate at column 66 for PDB format
                    out_file.write(pdb_line)
        
        # Copy ligand atoms with unique residue IDs
        res_id = 1
        for ligand_file in ligand_pdbqt_files:
            with ligand_file.open('r', encoding='utf-8') as lig_file:
                for line in lig_file:
                    if line.startswith(('ATOM', 'HETATM')):
                        # Update residue ID to ensure uniqueness
                        modified_line = line[:17] + f"{res_id:4d}" + line[22:66] + '\n'
                        out_file.write(modified_line)
                res_id += 1
    
    print(f"Created merged complex file with {len(ligand_pdbqt_files)} ligands")


def main(argv: List[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default=str(CONFIG_PATH), help="Path to YAML config file")
    parser.add_argument("--dry-run", action="store_true", help="Print commands without executing")
    args = parser.parse_args(argv)

    config = load_config(Path(args.config))
    scripts_dir = Path(__file__).resolve().parent
    working_dir = (scripts_dir / config["paths"]["working_dir"]).resolve()
    template_dir = (working_dir / config["wrapper"]["template_dir"]).resolve()

    receptor_pdbqt = working_dir / config["inputs"]["receptor_pdbqt"]
    ligand_pdbqt = working_dir / config["inputs"]["ligand_pdbqt"]
    output_dir = working_dir / config["wrapper"]["output_dir"]
    output_dir.mkdir(parents=True, exist_ok=True)

    # Copy receptor and ligand pdbqt to output_dir for AutoGrid/AutoDock
    import shutil
    receptor_in_output = output_dir / receptor_pdbqt.name
    ligand_in_output = output_dir / ligand_pdbqt.name
    if receptor_pdbqt.exists():
        shutil.copy2(receptor_pdbqt, receptor_in_output)
    if ligand_pdbqt.exists():
        shutil.copy2(ligand_pdbqt, ligand_in_output)

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
        "receptor": receptor_pdbqt.name,
        "center_x": str(config["autogrid"]["center"][0]),
        "center_y": str(config["autogrid"]["center"][1]),
        "center_z": str(config["autogrid"]["center"][2]),
        "spacing": str(config["autogrid"]["spacing"]),
        "ligand_types": config["inputs"]["ligand_types"],
        "map_definitions": "\n".join(map_definitions),
    }

    gpf_content = render_template(template_dir / "gpf_template.txt", gpf_mapping)
    write_file(gpf_path, gpf_content)
    
    # Copy AD4_parameters.dat to output directory for AutoGrid
    param_file_src = template_dir.parent / "AD4_parameters.dat"
    if param_file_src.exists():
        param_file_dst = output_dir / "AD4_parameters.dat"
        shutil.copy2(param_file_src, param_file_dst)
        print(f"Copied AD4_parameters.dat to {param_file_dst}")
    else:
        print(f"[WARN] AD4_parameters.dat not found at {param_file_src}")

    autogrid_exe = config["paths"]["autogrid4"]
    # Determine if we need WSL: command not on Windows but exists in WSL
    use_wsl_autogrid = not windows_command_exists(autogrid_exe) and wsl_command_exists(autogrid_exe)
    
    # Use relative paths in the output directory
    autogrid_cmd = build_command(
        autogrid_exe,
        [
            "-p",
            gpf_path.name,
            "-l",
            "autogrid.log",
        ],
        use_wsl=use_wsl_autogrid,
    )
    # Run autogrid in output_dir (need to cd first for WSL)
    if use_wsl_autogrid:
        wsl_output_dir = to_wsl_path(output_dir)
        autogrid_cmd = ["wsl", "bash", "-c", f"cd {wsl_output_dir} && {autogrid_exe} -p {gpf_path.name} -l autogrid.log"]
    run_command(autogrid_cmd, args.dry_run, cwd=output_dir if not use_wsl_autogrid else None)

    dpf_template = template_dir / "dpf_template.txt"
    autodock_exe = config["paths"]["autodock4"]
    # Determine if we need WSL for autodock
    use_wsl_autodock = not windows_command_exists(autodock_exe) and wsl_command_exists(autodock_exe)

    for seed in config["wrapper"]["seeds"]:
        map_definitions = "\n".join([
            f"map {ligand_type}.map" for ligand_type in config["inputs"]["ligand_types"].split()
        ] + ["elecmap e.map", "dsolvmap d.map"])
        mapping = {
            "seed": str(seed),
            "ligand_types": config["inputs"]["ligand_types"],
            "gridfld": gridfld.name,
            "maps": map_definitions,
            "ligand_pdbqt": ligand_pdbqt.name,
            "receptor": receptor_pdbqt.stem,
            "center_x": str(config["autogrid"]["center"][0]),
            "center_y": str(config["autogrid"]["center"][1]),
            "center_z": str(config["autogrid"]["center"][2]),
        }
        dpf_content = render_template(dpf_template, mapping)
        dpf_path = output_dir / f"wrapper_{seed}.dpf"
        write_file(dpf_path, dpf_content)
        dlg_path = output_dir / f"wrapper_{seed}.dlg"

        # Run autodock in output_dir
        if use_wsl_autodock:
            wsl_output_dir = to_wsl_path(output_dir)
            cmd = ["wsl", "bash", "-c", f"cd {wsl_output_dir} && {autodock_exe} -p {dpf_path.name} -l {dlg_path.name}"]
        else:
            cmd = build_command(
                autodock_exe,
                ["-p", dpf_path.name, "-l", dlg_path.name],
                use_wsl=False,
            )
        run_command(cmd, args.dry_run, cwd=output_dir if not use_wsl_autodock else None)

    # After all docking runs, merge the results into a single complex file
    print("\n=== Merging docking results ===")
    
    # Find all generated DLG files and extract best poses
    dlg_files = list(output_dir.glob("wrapper_*.dlg"))
    ligand_pdbqt_files = []
    
    for dlg_file in dlg_files:
        # Extract ligand ID from filename
        ligand_id = dlg_file.stem.replace("wrapper_", "")
        output_ligand = output_dir / f"best_pose_{ligand_id}.pdbqt"
        
        # Extract best pose from DLG file
        try:
            extract_best_pose(dlg_file, output_ligand)
            ligand_pdbqt_files.append(output_ligand)
        except Exception as e:
            print(f"Warning: Could not extract pose from {dlg_file}: {e}")
    
    # Create merged complex file
    if ligand_pdbqt_files:
        wrapped_complex_pdb = output_dir / "wrapped_complex.pdb"
        merge_complex_files(receptor_in_output, ligand_pdbqt_files, wrapped_complex_pdb)
        print(f"Successfully created {wrapped_complex_pdb} with {len(ligand_pdbqt_files)} ligands")
    else:
        print("Warning: No ligand poses were extracted, skipping complex merge")


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


if __name__ == "__main__":
    main()
