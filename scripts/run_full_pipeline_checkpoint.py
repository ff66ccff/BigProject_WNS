#!/usr/bin/env python3
"""Automate the entire docking and MD workflow from raw inputs with checkpoint support."""

from __future__ import annotations

import argparse
import shlex
import subprocess
import sys
from pathlib import Path
from typing import Any, Dict, Iterable, List

try:
    import yaml  # type: ignore
except ImportError as exc:  # pragma: no cover
    raise SystemExit("PyYAML is required. Install it via 'pip install pyyaml'.") from exc

SCRIPTS_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPTS_DIR.parent
DEFAULT_CONFIG = SCRIPTS_DIR / "config.yml"
UPDATE_MANIFEST = SCRIPTS_DIR / "update_manifest.py"

# Import state management
sys.path.append(str(SCRIPTS_DIR.parent / 'utils'))
from state_manager import StateManager


class PipelineError(RuntimeError):
    """Domain-specific error for pipeline failures."""


def load_config(path: Path) -> Dict[str, Any]:
    if not path.exists():
        raise PipelineError(
            f"Configuration file '{path}' not found. Copy config.example.yml and edit it first."
        )
    with path.open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def ensure_exists(path: Path, description: str) -> None:
    if not path.exists():
        raise PipelineError(f"Expected {description} at '{path}' but it was not found.")


def make_parent(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)


def run_command(label: str, cmd: List[str], dry_run: bool, cwd: Path | None = None) -> bool:
    prefix = f"[STEP] {label}:"
    print(prefix, " ".join(shlex.quote(part) for part in cmd))
    if dry_run:
        return False
    result = subprocess.run(cmd, cwd=cwd, check=False)
    if result.returncode != 0:
        raise PipelineError(f"Command for '{label}' failed with exit code {result.returncode}.")
    return True


def run_wsl(label: str, command: str, dry_run: bool) -> bool:
    cmd = ["wsl", "bash", "-lc", command]
    return run_command(label, cmd, dry_run)


def wsl_command_exists(cmd_name: str, conda_env: str | None = None) -> bool:
    """Check whether a command exists in WSL PATH by invoking `which`."""
    try:
        if conda_env:
            check_cmd = f"source ~/miniconda3/etc/profile.d/conda.sh && conda activate {conda_env} && which {cmd_name}"
            result = subprocess.run(["wsl", "bash", "-c", check_cmd], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        else:
            result = subprocess.run(["wsl", "which", cmd_name], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        return result.returncode == 0
    except Exception:
        return False


def to_wsl_path(path: Path) -> str:
    path = path.resolve()
    if path.drive:
        drive_letter = path.drive.rstrip(":").lower()
        return f"/mnt/{drive_letter}{path.as_posix()[2:]}"
    return path.as_posix()


def rel_to_workdir(path: Path, workdir: Path) -> str:
    return path.resolve().relative_to(workdir.resolve()).as_posix()


def update_manifest(
    python_exe: str,
    manifest_path: Path,
    commands: Iterable[str] | None = None,
    seeds: Iterable[str] | None = None,
    inputs: Dict[str, str] | None = None,
    outputs: Dict[str, str] | None = None,
    notes: str | None = None,
) -> None:
    args: List[str] = [python_exe, str(UPDATE_MANIFEST), "--manifest", str(manifest_path)]
    for command in commands or []:
        args.extend(["--command", command])
    for seed in seeds or []:
        args.extend(["--seed", seed])
    for key, value in (inputs or {}).items():
        args.extend(["--input", f"{key}={value}"])
    for key, value in (outputs or {}).items():
        args.extend(["--output", f"{key}={value}"])
    if notes is not None:
        args.extend(["--append-note", notes])
    subprocess.run(args, check=True)


def build_mgl_paths(mgl_root: Path) -> Dict[str, Path]:
    candidate_bases = [
        mgl_root / "MGLTools-1.5.7",
        mgl_root / "mgltools_x86_64_1.5.7",
        mgl_root,
    ]

    for base in candidate_bases:
        python_path = base / "mglpython.exe"
        if not python_path.exists():
            python_path = base / "python.exe"
        receptor_script = (
            base
            / "Lib"
            / "site-packages"
            / "AutoDockTools"
            / "Utilities24"
            / "prepare_receptor4.py"
        )
        ligand_script = (
            base
            / "Lib"
            / "site-packages"
            / "AutoDockTools"
            / "Utilities24"
            / "prepare_ligand4.py"
        )
        if python_path.exists() and receptor_script.exists() and ligand_script.exists():
            return {
                "python": python_path,
                "prepare_receptor": receptor_script,
                "prepare_ligand": ligand_script,
            }

    raise PipelineError(
        "Unable to locate a valid MGLTools installation. Verify 'paths.mgltools_root' in config.yml."
    )


def resolve_pose_files(working_dir: Path, entries: Iterable[str]) -> List[Path]:
    pose_paths: List[Path] = []
    for entry in entries:
        pattern = entry.replace("\\", "/")
        if any(char in pattern for char in "*?[]"):
            matches = list(working_dir.glob(pattern))
            pose_paths.extend(sorted(matches))
        else:
            pose_paths.append((working_dir / pattern).resolve())
    unique = []
    seen = set()
    for path in pose_paths:
        if path.exists() and path not in seen:
            unique.append(path)
            seen.add(path)
    return unique


def main(argv: List[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default=str(DEFAULT_CONFIG), help="Path to config.yml")
    parser.add_argument("--dry-run", action="store_true", help="Print commands without executing")
    parser.add_argument("--reset", action="store_true", help="Reset all progress and start from scratch")
    args = parser.parse_args(argv)

    config = load_config(Path(args.config))
    working_dir = (SCRIPTS_DIR / config["paths"]["working_dir"]).resolve()
    python_exe = config["paths"].get("python", "python")
    manifest_path = (working_dir / config.get("manifest", {}).get("path", "manifest/run-manifest.yml")).resolve()
    
    # Initialize state manager
    state_file = working_dir / "pipeline_state.json"
    state = StateManager(str(state_file))
    
    # Reset state if requested
    if args.reset:
        print("Resetting pipeline state...")
        state.data = {
            "step": 0,
            "current_stage": None,
            "completed_stages": []
        }
        state.update("reset", True)

    mgl_root = Path(config["paths"]["mgltools_root"]).resolve()
    mgl_paths = build_mgl_paths(mgl_root)
    ensure_exists(mgl_paths["python"], "mglpython.exe")
    ensure_exists(mgl_paths["prepare_receptor"], "prepare_receptor4.py")
    ensure_exists(mgl_paths["prepare_ligand"], "prepare_ligand4.py")

    raw_protein = (working_dir / config["inputs"]["raw_protein_pdb"]).resolve()
    cleaned_protein = (working_dir / config["inputs"]["cleaned_protein_pdb"]).resolve()
    receptor_pdbqt = (working_dir / config["inputs"]["receptor_pdbqt"]).resolve()
    raw_ligand = (working_dir / config["inputs"]["raw_ligand_mol2"]).resolve()
    ligand_pdbqt = (working_dir / config["inputs"]["ligand_pdbqt"]).resolve()
    ligand_basename = (working_dir / config["inputs"]["ligand_basename"]).resolve()
    complex_output = (working_dir / config["inputs"]["complex_output_pdb"]).resolve()

    ensure_exists(raw_protein, "raw protein PDB input")
    ensure_exists(raw_ligand, "raw ligand input")

    receptor_seed = str(config.get("receptor", {}).get("seed", "20231129"))
    wrapper_seeds = [str(seed) for seed in config.get("wrapper", {}).get("seeds", [])]
    
    # Define pipeline stages
    pipeline_stages = [
        "clean_protein",
        "prepare_receptor",
        "prepare_ligand", 
        "parameterize_ligand",
        "run_wrapper",
        "build_complex",
        "run_gromacs"
    ]
    
    # Get current stage
    current_stage = state.get("current_stage", None)
    completed_stages = state.get("completed_stages", [])
    
    print(f"Current stage: {current_stage}")
    print(f"Completed stages: {completed_stages}")
    
    # Determine starting stage
    if current_stage and current_stage in pipeline_stages:
        start_idx = pipeline_stages.index(current_stage)
    else:
        start_idx = 0
    
    # Run stages from starting point
    for i, stage in enumerate(pipeline_stages[start_idx:], start=start_idx):
        if stage in completed_stages:
            print(f"Skipping completed stage: {stage}")
            continue
            
        print(f"\n=== Running stage: {stage} ===")
        
        try:
            if stage == "clean_protein":
                # Step 1: clean protein PDB
                make_parent(cleaned_protein)
                preprocess_cmd = [
                    python_exe,
                    str(SCRIPTS_DIR / "preprocess_pdb.py"),
                    str(raw_protein),
                    "-o",
                    str(cleaned_protein),
                ]
                executed = run_command("Clean protein PDB", preprocess_cmd, args.dry_run)
                if executed:
                    update_manifest(
                        python_exe,
                        manifest_path,
                        commands=[" ".join(preprocess_cmd)],
                        inputs={"raw_protein": str(raw_protein)},
                        outputs={"cleaned_protein": str(cleaned_protein)},
                    )
                state.update("current_stage", "prepare_receptor")
                
            elif stage == "prepare_receptor":
                # Step 2: receptor PDBQT via MGLTools
                make_parent(receptor_pdbqt)
                receptor_cmd = [
                    str(mgl_paths["python"]),
                    str(mgl_paths["prepare_receptor"]),
                    "-r",
                    cleaned_protein.name,
                    "-o",
                    str(receptor_pdbqt),
                    "-A",
                    "checkhydrogens",
                    "-U",
                    "nphs_lps",
                    "-e",
                    f"seed={receptor_seed}",
                ]
                executed = run_command(
                    "Prepare receptor PDBQT",
                    receptor_cmd,
                    args.dry_run,
                    cwd=cleaned_protein.parent,
                )
                if executed:
                    update_manifest(
                        python_exe,
                        manifest_path,
                        commands=[" ".join(receptor_cmd)],
                        seeds=[receptor_seed],
                        outputs={"receptor_pdbqt": str(receptor_pdbqt)},
                    )
                state.update("current_stage", "prepare_ligand")
                
            elif stage == "prepare_ligand":
                # Step 3: ligand PDBQT via MGLTools
                make_parent(ligand_pdbqt)
                ligand_cmd = [
                    str(mgl_paths["python"]),
                    str(mgl_paths["prepare_ligand"]),
                    "-l",
                    raw_ligand.name,
                    "-o",
                    str(ligand_pdbqt),
                    "-A",
                    "checkhydrogens",
                ]
                executed = run_command(
                    "Prepare ligand PDBQT",
                    ligand_cmd,
                    args.dry_run,
                    cwd=raw_ligand.parent,
                )
                if executed:
                    update_manifest(
                        python_exe,
                        manifest_path,
                        commands=[" ".join(ligand_cmd)],
                        outputs={"ligand_pdbqt": str(ligand_pdbqt)},
                    )
                state.update("current_stage", "parameterize_ligand")
                
            elif stage == "parameterize_ligand":
                # Step 4: AmberTools / ACPYPE parameterization (WSL)
                ligand_script = (working_dir / config["ambertools"]["ligand_script"]).resolve()
                ensure_exists(ligand_script, "ligand_param.sh script")
                make_parent(ligand_basename)
                wsl_root = to_wsl_path(working_dir)

                amber_cfg = config.get("ambertools", {})
                skip_amber = bool(amber_cfg.get("skip", False))
                fail_on_missing = bool(amber_cfg.get("fail_on_missing", False))

                ligand_param_cmd = (
                    f"cd {shlex.quote(wsl_root)} && "
                    f"bash {shlex.quote(rel_to_workdir(ligand_script, working_dir))} "
                    f"{shlex.quote(rel_to_workdir(raw_ligand, working_dir))} "
                    f"{shlex.quote(rel_to_workdir(ligand_basename, working_dir))}"
                )

                if skip_amber:
                    print("[INFO] Skipping AmberTools ligand parameterization (config: ambertools.skip = true)")
                    update_manifest(
                        python_exe,
                        manifest_path,
                        notes="Skipped AmberTools ligand parameterization (config: ambertools.skip = true)",
                    )
                else:
                    executed = run_wsl("AmberTools ligand parameterization", ligand_param_cmd, args.dry_run)
                    if executed:
                        update_manifest(
                            python_exe,
                            manifest_path,
                            commands=[ligand_param_cmd],
                            inputs={"raw_ligand": str(raw_ligand)},
                            outputs={"ligand_basename": str(ligand_basename)},
                        )
                state.update("current_stage", "run_wrapper")
                
            elif stage == "run_wrapper":
                # Step 5: Wrapper docking
                wrapper_cfg = config.get("wrapper", {})
                wrapper_script = (working_dir / wrapper_cfg["script"]).resolve()
                ensure_exists(wrapper_script, "wrapper script")
                wrapper_output_dir = (working_dir / wrapper_cfg["output_dir"]).resolve()
                make_parent(wrapper_output_dir / "placeholder")

                wrapper_cmd = [
                    python_exe,
                    str(wrapper_script),
                    "--receptor",
                    str(receptor_pdbqt),
                    "--ligand",
                    str(ligand_pdbqt),
                    "--output-dir",
                    str(wrapper_output_dir),
                ]

                for seed in wrapper_seeds:
                    wrapper_cmd.extend(["--seed", seed])

                if "center" in wrapper_cfg:
                    wrapper_cmd.extend(["--center"] + [str(x) for x in wrapper_cfg["center"]])
                if "size" in wrapper_cfg:
                    wrapper_cmd.extend(["--size"] + [str(x) for x in wrapper_cfg["size"]])

                executed = run_command("Wrapper docking", wrapper_cmd, args.dry_run)
                if executed:
                    update_manifest(
                        python_exe,
                        manifest_path,
                        commands=[" ".join(wrapper_cmd)],
                        inputs={"receptor_pdbqt": str(receptor_pdbqt), "ligand_pdbqt": str(ligand_pdbqt)},
                        outputs={"wrapper_output_dir": str(wrapper_output_dir)},
                        seeds=wrapper_seeds,
                    )
                state.update("current_stage", "build_complex")
                
            elif stage == "build_complex":
                # Step 6: Build complex
                wrapper_cfg = config.get("wrapper", {})
                wrapper_output_dir = (working_dir / wrapper_cfg["output_dir"]).resolve()
                pose_entries = wrapper_cfg.get("poses", ["docked_ligand_*.pdbqt"])
                pose_paths = resolve_pose_files(wrapper_output_dir, pose_entries)

                if not pose_paths:
                    if wrapper_cfg.get("fail_on_missing", True):
                        raise PipelineError("No docked poses found after wrapper stage")
                    else:
                        print("[WARNING] No docked poses found, skipping complex building")
                        state.update("current_stage", "run_gromacs")
                        continue

                make_parent(complex_output)
                complex_cmd = [
                    python_exe,
                    str(SCRIPTS_DIR / "build_complex.py"),
                    str(cleaned_protein),
                    *[str(path) for path in pose_paths],
                    "-o",
                    str(complex_output),
                ]
                executed = run_command("Build complex", complex_cmd, args.dry_run)
                if executed:
                    update_manifest(
                        python_exe,
                        manifest_path,
                        commands=[" ".join(complex_cmd)],
                        outputs={"complex_pdb": str(complex_output)},
                    )
                state.update("current_stage", "run_gromacs")
                
            elif stage == "run_gromacs":
                # Step 7: GROMACS pipeline (WSL) - 使用新的 Python 脚本
                gromacs_script = (SCRIPTS_DIR / "gromacs_pipeline.py").resolve()
                ensure_exists(gromacs_script, "gromacs pipeline script")
                
                gmx_workdir = (working_dir / config["gromacs"]["workdir"]).resolve()
                make_parent(gmx_workdir / "placeholder")
                
                # 使用 Python 脚本而不是 bash 脚本
                gromacs_cmd = [
                    python_exe,
                    str(gromacs_script),
                    str(complex_output),
                    str(gmx_workdir),
                    str(raw_ligand)
                ]
                
                executed = run_command("GROMACS pipeline", gromacs_cmd, args.dry_run)
                if executed:
                    update_manifest(
                        python_exe,
                        manifest_path,
                        commands=[" ".join(gromacs_cmd)],
                        outputs={"gromacs_workdir": str(gmx_workdir)},
                    )
                state.update("current_stage", "completed")
            
            # Mark stage as completed
            if stage not in completed_stages:
                completed_stages.append(stage)
                state.update("completed_stages", completed_stages)
                
        except Exception as e:
            print(f"Error in stage {stage}: {e}")
            state.update("current_stage", stage)  # Save current stage for retry
            raise

    print("\n=== Pipeline completed successfully! ===")


if __name__ == "__main__":
    try:
        main(sys.argv[1:])
    except PipelineError as exc:
        raise SystemExit(str(exc))