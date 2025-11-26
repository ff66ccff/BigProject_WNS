#!/usr/bin/env python3
"""Create a reproducibility archive with manifests, scripts, and outputs."""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path
from typing import Iterable


DEFAULT_INCLUDE = [
    "manifest/run-manifest.yml",
    "scripts",
    "analysis",
    "autodock_runs",
    "gmx",
    "md",
]


def collect_files(root: Path, include: Iterable[str], output: Path) -> None:
    if output.exists():
        output.unlink()
    temp_dir = root / "_repro_pack"
    if temp_dir.exists():
        shutil.rmtree(temp_dir)
    temp_dir.mkdir()
    try:
        for item in include:
            source = root / item
            if source.exists():
                if source.is_dir():
                    shutil.copytree(source, temp_dir / source.name)
                else:
                    shutil.copy2(source, temp_dir / source.name)
        shutil.make_archive(str(output), "zip", root_dir=str(temp_dir))
    finally:
        shutil.rmtree(temp_dir, ignore_errors=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(".."), help="Project root directory")
    parser.add_argument("--output", type=Path, default=Path("../manifest/repro_pack"), help="Archive output base path")
    args = parser.parse_args()

    collect_files(args.root.resolve(), DEFAULT_INCLUDE, args.output.resolve())


if __name__ == "__main__":
    main()
