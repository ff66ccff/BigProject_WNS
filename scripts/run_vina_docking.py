#!/usr/bin/env python3
"""Simple Vina docking for WnS blind surface docking."""

from __future__ import annotations

import argparse
import subprocess
from pathlib import Path

try:
    import yaml
except ImportError:
    raise SystemExit("PyYAML required: pip install pyyaml")


def to_wsl_path(path: Path) -> str:
    """Convert Windows path to WSL path."""
    path_str = str(path)
    if ':' in path_str:
        import os
        drive, tail = os.path.splitdrive(path_str)
        if drive:
            drive_letter = drive[0].lower()
            linux_path = tail.replace('\\', '/').lstrip('/')
            return f"/mnt/{drive_letter}/{linux_path}"
    return path_str.replace('\\', '/')


def run_vina_docking(config: dict) -> None:
    """Run Vina docking."""
    scripts_dir = Path(__file__).resolve().parent
    working_dir = (scripts_dir / config["paths"]["working_dir"]).resolve()
    
    receptor = working_dir / config["inputs"]["receptor_pdbqt"]
    ligand = working_dir / config["inputs"]["ligand_pdbqt"]
    output_dir = working_dir / config["wrapper"]["output_dir"]
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Grid parameters
    center = config["autogrid"]["center"]
    # Calculate box size from npts and spacing
    npts = config["autogrid"]["npts"]
    spacing = config["autogrid"]["spacing"]
    size = [n * spacing for n in npts]
    
    seeds = config["wrapper"]["seeds"]
    
    for seed in seeds:
        output_pdbqt = output_dir / f"vina_docked_{seed}.pdbqt"
        log_file = output_dir / f"vina_{seed}.log"
        
        print(f"\n=== Vina docking (seed: {seed}) ===")
        
        cmd = [
            "wsl", "vina",
            "--receptor", to_wsl_path(receptor),
            "--ligand", to_wsl_path(ligand),
            "--center_x", str(center[0]),
            "--center_y", str(center[1]),
            "--center_z", str(center[2]),
            "--size_x", str(size[0]),
            "--size_y", str(size[1]),
            "--size_z", str(size[2]),
            "--seed", str(seed),
            "--exhaustiveness", "8",
            "--num_modes", "9",
            "--out", to_wsl_path(output_pdbqt),
        ]
        
        print(" ".join(cmd))
        result = subprocess.run(cmd, check=False)
        
        if result.returncode == 0:
            print(f"Success! Output: {output_pdbqt}")
        else:
            print(f"Failed with exit code {result.returncode}")
            break
    
    print("\n=== Vina docking completed ===")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="config.yml", help="YAML config file")
    args = parser.parse_args()
    
    with open(args.config, "r", encoding="utf-8") as f:
        config = yaml.safe_load(f)
    
    run_vina_docking(config)


if __name__ == "__main__":
    main()
