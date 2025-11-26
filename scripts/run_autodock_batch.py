#!/usr/bin/env python3
"""Generate AutoGrid/AutoDock input files and run batch docking jobs."""

from __future__ import annotations

import argparse
import os
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


def build_command(executable: str, args: List[str]) -> List[str]:
    if executable.startswith("/"):
        return ["wsl", executable, *args]
    return [executable, *args]


def run_command(cmd: List[str], dry_run: bool) -> None:
    print(" ".join(cmd))
    if dry_run:
        return
    result = subprocess.run(cmd, check=False)
    if result.returncode != 0:
        raise RuntimeError(f"Command '{cmd}' failed with exit code {result.returncode}")


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

    gpf_path = output_dir / "autogrid.gpf"
    gridfld = output_dir / "protein.maps.fld"

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
    }

    gpf_content = render_template(template_dir / "gpf_template.txt", gpf_mapping)
    write_file(gpf_path, gpf_content)

    autogrid_exe = config["paths"]["autogrid4"]
    autogrid_cmd = build_command(
        autogrid_exe,
        [
            "-p",
            to_wsl_path(gpf_path) if autogrid_exe.startswith("/") else str(gpf_path),
            "-l",
            to_wsl_path(output_dir / "autogrid.log")
            if autogrid_exe.startswith("/")
            else str(output_dir / "autogrid.log"),
        ],
    )
    run_command(autogrid_cmd, args.dry_run)

    dpf_template = template_dir / "dpf_template.txt"
    autodock_exe = config["paths"]["autodock4"]

    for seed in config["wrapper"]["seeds"]:
        maps = " ".join([
            f"{ligand_type}.map" for ligand_type in config["inputs"]["ligand_types"].split()
        ])
        mapping = {
            "seed": str(seed),
            "ligand_types": config["inputs"]["ligand_types"],
            "gridfld": gridfld.name,
            "maps": maps,
            "ligand_pdbqt": ligand_pdbqt.name,
            "center_x": str(config["autogrid"]["center"][0]),
            "center_y": str(config["autogrid"]["center"][1]),
            "center_z": str(config["autogrid"]["center"][2]),
        }
        dpf_content = render_template(dpf_template, mapping)
        dpf_path = output_dir / f"wrapper_{seed}.dpf"
        write_file(dpf_path, dpf_content)
        dlg_path = output_dir / f"wrapper_{seed}.dlg"

        cmd = build_command(
            autodock_exe,
            [
                "-p",
                to_wsl_path(dpf_path) if autodock_exe.startswith("/") else str(dpf_path),
                "-l",
                to_wsl_path(dlg_path) if autodock_exe.startswith("/") else str(dlg_path),
            ],
        )
        run_command(cmd, args.dry_run)


if __name__ == "__main__":
    main()
