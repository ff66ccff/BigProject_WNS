#!/usr/bin/env python3
"""Wrap 'n' Shake docking pipeline: sequential docking with atom masking.

Supports checkpointing - can be interrupted and resumed from the last completed iteration.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import shutil
import subprocess
import sys
from pathlib import Path
from string import Template
from typing import Dict, List, Optional

try:
    import yaml  # type: ignore
except ImportError as exc:  # pragma: no cover
    raise SystemExit("PyYAML is required. Install it via 'pip install pyyaml'.") from exc

# Import our masking function
from mask_pdbqt import mask_receptor


class CheckpointState:
    """Manages checkpoint state for resumable docking."""
    
    def __init__(self, checkpoint_file: Path):
        self.checkpoint_file = checkpoint_file
        self.state = {
            "completed_seeds": [],
            "docked_ligand_files": [],
            "current_receptor": None,
            "autogrid_complete": False,
            "successful_docks": 0,
        }
    
    def load(self) -> bool:
        """Load checkpoint state from file. Returns True if loaded successfully."""
        if self.checkpoint_file.exists():
            try:
                with self.checkpoint_file.open("r", encoding="utf-8") as f:
                    self.state = json.load(f)
                print(f"[CHECKPOINT] Resuming from checkpoint: {len(self.state['completed_seeds'])} seeds completed")
                return True
            except (json.JSONDecodeError, KeyError) as e:
                print(f"[CHECKPOINT] Warning: Could not load checkpoint ({e}), starting fresh")
                return False
        return False
    
    def save(self) -> None:
        """Save current state to checkpoint file."""
        with self.checkpoint_file.open("w", encoding="utf-8") as f:
            json.dump(self.state, f, indent=2)
        print(f"[CHECKPOINT] State saved: {len(self.state['completed_seeds'])} seeds completed")
    
    def mark_autogrid_complete(self) -> None:
        self.state["autogrid_complete"] = True
        self.save()
    
    def is_autogrid_complete(self) -> bool:
        return self.state.get("autogrid_complete", False)
    
    def is_seed_completed(self, seed: int) -> bool:
        return seed in self.state["completed_seeds"]
    
    def mark_seed_completed(self, seed: int, ligand_file: Optional[str], receptor_file: str) -> None:
        self.state["completed_seeds"].append(seed)
        if ligand_file:
            self.state["docked_ligand_files"].append(ligand_file)
            self.state["successful_docks"] = len(self.state["docked_ligand_files"])
        self.state["current_receptor"] = receptor_file
        self.save()
    
    def get_current_receptor(self) -> Optional[str]:
        return self.state.get("current_receptor")
    
    def get_docked_ligands(self) -> List[str]:
        return self.state.get("docked_ligand_files", [])
    
    def get_successful_docks(self) -> int:
        return self.state.get("successful_docks", 0)
    
    def reset(self) -> None:
        """Reset checkpoint state."""
        self.state = {
            "completed_seeds": [],
            "docked_ligand_files": [],
            "current_receptor": None,
            "autogrid_complete": False,
            "successful_docks": 0,
        }
        if self.checkpoint_file.exists():
            self.checkpoint_file.unlink()
        print("[CHECKPOINT] State reset")


def to_wsl_path(path: Path) -> str:
    """Convert Windows path to WSL path."""
    path_str = str(path)
    
    # Handle absolute Windows paths (e.g., D:\Project\data)
    if ':' in path_str:
        import os
        drive, tail = os.path.splitdrive(path_str)
        if drive:
            drive_letter = drive[0].lower()
            # Convert backslashes to forward slashes and remove leading backslash
            linux_path = tail.replace('\\', '/').lstrip('/')
            return f"/mnt/{drive_letter}/{linux_path}"
    
    # Handle relative paths or already Unix-style paths
    return path_str.replace('\\', '/')


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
        if "DOCKED: USER    Run = 1" in line:
            # Found best pose, extract coordinates
            j = i
            while j < len(lines) and not lines[j].startswith("DOCKED: USER    Run = 2"):
                if lines[j].startswith("DOCKED: ATOM"):
                    best_pose_lines.append(lines[j][8:])  # Remove "DOCKED: " prefix
                j += 1
            found_best = True
            break
    
    if not found_best:
        raise RuntimeError(f"Could not find best pose in {dlg_file}")
    
# Write best pose to PDBQT file
    with output_pdbqt.open('w', encoding='utf-8') as f:
        # Write PDBQT header
        f.write("TORSDOF 0\n")
        # Write atoms
        for line in best_pose_lines:
            f.write(line)
    
    print(f"Extracted best pose to {output_pdbqt}")


def check_ligand_clash(new_ligand_path: Path, existing_ligands: List[Path], 
                       min_distance: float = 2.0) -> bool:
    """Check if new ligand clashes with existing ligands."""
    if not existing_ligands:
        return False
    
    # Load new ligand atoms
    new_atoms = []
    for line in new_ligand_path.read_text(encoding="utf-8").splitlines():
        if line.startswith(("ATOM", "HETATM")):
            try:
                coord = (
                    float(line[30:38]),
                    float(line[38:46]),
                    float(line[46:54]),
                )
                new_atoms.append(coord)
            except (ValueError, IndexError):
                continue
    
    if not new_atoms:
        return True  # Treat as clash if no valid atoms
    
    # Check against existing ligands
    for existing_ligand in existing_ligands:
        existing_atoms = []
        for line in existing_ligand.read_text(encoding="utf-8").splitlines():
            if line.startswith(("ATOM", "HETATM")):
                try:
                    coord = (
                        float(line[30:38]),
                        float(line[38:46]),
                        float(line[46:54]),
                    )
                    existing_atoms.append(coord)
                except (ValueError, IndexError):
                    continue
        
        # Calculate minimum distance
        for new_coord in new_atoms:
            for existing_coord in existing_atoms:
                dist = math.sqrt(
                    (new_coord[0] - existing_coord[0])**2 +
                    (new_coord[1] - existing_coord[1])**2 +
                    (new_coord[2] - existing_coord[2])**2
                )
                if dist < min_distance:
                    return True  # Clash detected
    
    return False  # No clash


def run_wrap_n_shake_docking(config: Dict, dry_run: bool = False, reset_checkpoint: bool = False) -> None:
    """Run Wrap 'n' Shake docking pipeline with checkpoint support.
    
    Args:
        config: Configuration dictionary
        dry_run: If True, only print commands without executing
        reset_checkpoint: If True, reset checkpoint and start fresh
    """
    scripts_dir = Path(__file__).resolve().parent
    working_dir = (scripts_dir / config["paths"]["working_dir"]).resolve()
    template_dir = (working_dir / config["wrapper"]["template_dir"]).resolve()
    
    receptor_pdbqt = working_dir / config["inputs"]["receptor_pdbqt"]
    ligand_pdbqt = working_dir / config["inputs"]["ligand_pdbqt"]
    output_dir = working_dir / config["wrapper"]["output_dir"]
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Initialize checkpoint system
    checkpoint_file = output_dir / "docking_checkpoint.json"
    checkpoint = CheckpointState(checkpoint_file)
    
    if reset_checkpoint:
        checkpoint.reset()
    else:
        checkpoint.load()
    
    # Determine current receptor (from checkpoint or original)
    if checkpoint.get_current_receptor():
        receptor_current = Path(checkpoint.get_current_receptor())
        print(f"[RESUME] Using masked receptor from checkpoint: {receptor_current.name}")
    else:
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
    
    # Run AutoGrid (skip if already complete from checkpoint)
    if checkpoint.is_autogrid_complete():
        print("[CHECKPOINT] AutoGrid already completed, skipping...")
    else:
        autogrid_exe = config["paths"]["autogrid4"]
        use_wsl_autogrid = not windows_command_exists(autogrid_exe) and wsl_command_exists(autogrid_exe)
        
        if use_wsl_autogrid:
            wsl_output_dir = to_wsl_path(output_dir)
            autogrid_cmd = ["wsl", "bash", "-c", f"cd {wsl_output_dir} && {autogrid_exe} -p {gpf_path.name} -l autogrid.log"]
        else:
            autogrid_cmd = build_command(autogrid_exe, ["-p", gpf_path.name, "-l", "autogrid.log"], use_wsl=False)
        
        run_command(autogrid_cmd, dry_run, cwd=output_dir if not use_wsl_autogrid else None)
        
        if not dry_run:
            checkpoint.mark_autogrid_complete()
    
    # Sequential docking loop with controls
    seeds = config["wrapper"]["seeds"]
    max_cycles = config.get("wrapper", {}).get("max_cycles", 20)
    min_ligand_distance = config.get("wrapper", {}).get("min_ligand_distance", 2.0)
    
    # Restore docked ligands from checkpoint
    docked_ligands = [Path(p) for p in checkpoint.get_docked_ligands()]
    successful_docks = checkpoint.get_successful_docks()
    
    for i, seed in enumerate(seeds):
        if successful_docks >= max_cycles:
            print(f"\n=== Maximum cycles ({max_cycles}) reached, stopping ===")
            break
        
        # Skip completed seeds from checkpoint
        if checkpoint.is_seed_completed(seed):
            print(f"\n=== Skipping seed {seed} (already completed from checkpoint) ===")
            continue
        
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
        
        # Skip post-processing in dry-run mode
        if dry_run:
            print(f"[DRY-RUN] Would extract best pose from {dlg_path}")
            continue
        
        # Extract best pose
        extract_best_pose(dlg_path, docked_ligand_path)
        
        # Check for ligand clashes
        if check_ligand_clash(docked_ligand_path, docked_ligands, min_ligand_distance):
            print(f"WARNING: Ligand {seed} clashes with existing ligands, discarding...")
            docked_ligand_path.unlink()  # Remove the clashed ligand
            # Save checkpoint even for failed ligands
            checkpoint.mark_seed_completed(seed, None, str(receptor_current))
            continue
        
        # Accept this ligand
        docked_ligands.append(docked_ligand_path)
        successful_docks += 1
        
        # Mask receptor atoms for next iteration
        if i < len(seeds) - 1 and successful_docks < max_cycles:
            print(f"Masking receptor atoms near docked ligand {seed}...")
            receptor_next = output_dir / f"receptor_masked_{seed}.pdbqt"
            mask_receptor(receptor_current, [docked_ligand_path], receptor_next, cutoff=3.5)
            receptor_current = receptor_next
        
        # Save checkpoint after successful docking
        checkpoint.mark_seed_completed(seed, str(docked_ligand_path), str(receptor_current))
    
    print(f"\n=== Wrap 'n' Shake docking completed ===")
    print(f"Successfully docked {len(docked_ligands)} ligands out of {len(seeds)} attempts:")
    for ligand_path in docked_ligands:
        print(f"  - {ligand_path}")


def main(argv: List[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="config.yml", help="Path to YAML config file")
    parser.add_argument("--dry-run", action="store_true", help="Print commands without executing")
    parser.add_argument("--reset", action="store_true", help="Reset checkpoint and start fresh")
    args = parser.parse_args(argv)
    
    config = load_config(Path(args.config))
    run_wrap_n_shake_docking(config, args.dry_run, args.reset)


if __name__ == "__main__":
    main()