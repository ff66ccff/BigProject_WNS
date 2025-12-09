#!/usr/bin/env python3
"""Automate entire docking and MD workflow from raw inputs."""

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
try:
    from state_manager import StateManager
except ImportError:
    # Fallback mock if utils not found/setup yet
    class StateManager:
        def __init__(self, path): self.data = {}
        def get(self, k, d=None): return self.data.get(k, d)
        def update(self, k, v): self.data[k] = v

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
    print(prefix, " ".join(shlex.quote(str(part)) for part in cmd))
    if dry_run:
        return False
    result = subprocess.run([str(c) for c in cmd], cwd=cwd, check=False)
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
            path = working_dir / pattern
            if path.exists():
                pose_paths.append(path.resolve())
    
    unique = []
    seen = set()
    for path in pose_paths:
        if path not in seen:
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
    make_parent(state_file)
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
    # Ensure MGL tools exist
    if not args.dry_run:
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

    if not args.dry_run:
        ensure_exists(raw_protein, "raw protein PDB input")
        ensure_exists(raw_ligand, "raw ligand input")

    receptor_seed = str(config.get("receptor", {}).get("seed", "20231129"))
    wrapper_seeds = [str(seed) for seed in config.get("wrapper", {}).get("seeds", [])]
    
    # Define pipeline stages in order
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
    completed_stages = state.get("completed_stages", []) or []
    
    print(f"Current recorded stage: {current_stage}")
    print(f"Completed stages: {completed_stages}")
    
    # Run stages
    for stage in pipeline_stages:
        if stage in completed_stages:
            print(f"Skipping completed stage: {stage}")
            continue
            
        print(f"\n=== Running stage: {stage} ===")
        
        try:
            # ----------------------------------------------------------------
            # STAGE: clean_protein
            # ----------------------------------------------------------------
            if stage == "clean_protein":
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
                        commands=[" ".join(map(str, preprocess_cmd))],
                        inputs={"raw_protein": str(raw_protein)},
                        outputs={"cleaned_protein": str(cleaned_protein)},
                    )

            # ----------------------------------------------------------------
            # STAGE: prepare_receptor
            # ----------------------------------------------------------------
            elif stage == "prepare_receptor":
                make_parent(receptor_pdbqt)
                receptor_cmd = [
                    str(mgl_paths["python"]),
                    str(mgl_paths["prepare_receptor"]),
                    "-r", cleaned_protein.name,
                    "-o", str(receptor_pdbqt),
                    "-A", "checkhydrogens",
                    "-U", "nphs_lps",
                    "-e", f"seed={receptor_seed}",
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
                        commands=[" ".join(map(str, receptor_cmd))],
                        seeds=[receptor_seed],
                        outputs={"receptor_pdbqt": str(receptor_pdbqt)},
                    )

            # ----------------------------------------------------------------
            # STAGE: prepare_ligand
            # ----------------------------------------------------------------
            elif stage == "prepare_ligand":
                make_parent(ligand_pdbqt)
                ligand_cmd = [
                    str(mgl_paths["python"]),
                    str(mgl_paths["prepare_ligand"]),
                    "-l", raw_ligand.name,
                    "-o", str(ligand_pdbqt),
                    "-A", "checkhydrogens",
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
                        commands=[" ".join(map(str, ligand_cmd))],
                        outputs={"ligand_pdbqt": str(ligand_pdbqt)},
                    )

            # ----------------------------------------------------------------
            # STAGE: parameterize_ligand (AmberTools)
            # ----------------------------------------------------------------
            elif stage == "parameterize_ligand":
                ligand_script = (working_dir / config["ambertools"]["ligand_script"]).resolve()
                if not args.dry_run:
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
                    print("[INFO] Skipping AmberTools parameterization (config: ambertools.skip=true)")
                else:
                    # Detect presence of antechamber in WSL
                    if not wsl_command_exists("antechamber", conda_env="amber"):
                        msg = "WSL 'antechamber' not found."
                        if fail_on_missing:
                            raise PipelineError(msg + " Install AmberTools or set skip=true.")
                        print(f"[WARN] {msg} Skipping parameterization.")
                    else:
                        executed = run_wsl("Ligand parameterization", ligand_param_cmd, args.dry_run)
                        if executed:
                            update_manifest(
                                python_exe,
                                manifest_path,
                                commands=[ligand_param_cmd],
                                outputs={"ligand_param_basename": str(ligand_basename)},
                            )

            # ----------------------------------------------------------------
            # STAGE: run_wrapper (AutoDock)
            # ----------------------------------------------------------------
            elif stage == "run_wrapper":
                autodock_cmd = [
                    python_exe,
                    str(SCRIPTS_DIR / "run_autodock_batch.py"),
                    "--config", str(Path(args.config).resolve()),
                ]
                if args.dry_run:
                    autodock_cmd.append("--dry-run")
                
                executed = run_command("AutoDock batch", autodock_cmd, args.dry_run)
                if executed:
                    update_manifest(
                        python_exe,
                        manifest_path,
                        commands=[" ".join(map(str, autodock_cmd))],
                        seeds=wrapper_seeds,
                    )

            # ----------------------------------------------------------------
            # STAGE: build_complex
            # ----------------------------------------------------------------
            elif stage == "build_complex":
                pose_entries = config["inputs"].get("pose_files", [])
                pose_paths = resolve_pose_files(working_dir, pose_entries)
                
                if not pose_paths:
                    print("[WARN] No pose files found. Skipping complex construction.")
                else:
                    make_parent(complex_output)
                    complex_cmd = [
                        python_exe,
                        str(SCRIPTS_DIR / "build_complex.py"),
                        str(cleaned_protein),
                        *[str(path) for path in pose_paths],
                        "-o", str(complex_output),
                    ]
                    executed = run_command("Build complex", complex_cmd, args.dry_run)
                    if executed:
                        update_manifest(
                            python_exe,
                            manifest_path,
                            commands=[" ".join(map(str, complex_cmd))],
                            outputs={"complex_pdb": str(complex_output)},
                        )

            # ----------------------------------------------------------------
            # STAGE: run_gromacs (Includes Filtering)
            # ----------------------------------------------------------------
            elif stage == "run_gromacs":
                gromacs_script = (working_dir / config["gromacs"].get("full_auto_script", "scripts/gromacs_full_auto.sh")).resolve()
                if not gromacs_script.exists():
                    gromacs_script = (working_dir / config["gromacs"]["pipeline_script"]).resolve()
                
                if not args.dry_run:
                    ensure_exists(gromacs_script, "gromacs script")
                
                gmx_workdir = (working_dir / config["gromacs"]["workdir"]).resolve()
                mdp_dir = (working_dir / config["gromacs"].get("mdp_dir", "mdp")).resolve()
                make_parent(gmx_workdir / "placeholder")
                wsl_root = to_wsl_path(working_dir)
                
                gmx_cmd = (
                    f"cd {shlex.quote(wsl_root)} && "
                    f"export GMX={shlex.quote(config['paths'].get('gmx', 'gmx'))} && "
                    f"bash {shlex.quote(rel_to_workdir(gromacs_script, working_dir))} "
                    f"{shlex.quote(rel_to_workdir(complex_output, working_dir))} "
                    f"{shlex.quote(rel_to_workdir(raw_ligand, working_dir))} "
                    f"{shlex.quote(rel_to_workdir(gmx_workdir, working_dir))} "
                    f"{shlex.quote(rel_to_workdir(mdp_dir, working_dir))}"
                )
                executed = run_wsl("GROMACS pipeline", gmx_cmd, args.dry_run)
                
                if executed:
                    update_manifest(
                        python_exe,
                        manifest_path,
                        commands=[gmx_cmd],
                        outputs={"gromacs_workdir": str(gmx_workdir)},
                    )
                    
                    # 7b. Post-Gromacs: Stable H-Bond Filtering
                    final_tpr = gmx_workdir / "final.tpr"
                    final_xtc = gmx_workdir / "final.xtc"
                    stable_hbond_output = working_dir / "final_results" / "stable_hbond_complex.pdb"
                    
                    if final_tpr.exists():
                        hbond_cmd = [
                            python_exe,
                            str(SCRIPTS_DIR / "filter_stable_hbonds.py"),
                            "--tpr", str(final_tpr),
                            "--output", str(stable_hbond_output)
                        ]
                        if final_xtc.exists():
                            hbond_cmd.extend(["--traj", str(final_xtc)])
                            
                        run_command("Stable H-Bond Filtering", hbond_cmd, args.dry_run)
                    else:
                        print("[WARN] final.tpr not found, skipping H-bond filter.")

            # Update state after success
            if not args.dry_run:
                completed_stages.append(stage)
                state.update("completed_stages", completed_stages)
                state.update("current_stage", stage)

        except Exception as e:
            print(f"❌ Error in stage '{stage}': {e}")
            raise PipelineError(f"Pipeline failed at stage {stage}") from e

    print("\n✅ Pipeline completed successfully.")


if __name__ == "__main__":
    try:
        main(sys.argv[1:])
    except PipelineError as exc:
        raise SystemExit(str(exc))