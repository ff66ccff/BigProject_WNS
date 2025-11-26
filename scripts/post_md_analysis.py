#!/usr/bin/env python3
"""Post-processing utilities for MD trajectories using MDAnalysis."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

try:
    import MDAnalysis as mda  # type: ignore
    from MDAnalysis.analysis import contacts, density, rms
except ImportError as exc:  # pragma: no cover
    raise SystemExit("Install MDAnalysis to run post_md_analysis.py") from exc


def occupancy(universe: mda.Universe, selection: str, output: Path) -> None:
    analysis = density.occupancy.Occupancy(universe, selection)
    analysis.run()
    analysis.save(str(output))


def contact_frequency(universe: mda.Universe, sel1: str, sel2: str, radius: float, output: Path) -> None:
    analysis = contacts.Contacts(
        universe,
        select=(sel1, sel2),
        refgroup=(sel1, sel2),
        radius=radius,
    )
    analysis.run()
    output.write_text(json.dumps(analysis.timeseries.tolist(), indent=2), encoding="utf-8")


def compute_rmsd(universe: mda.Universe, selection: str, reference: str, output: Path) -> None:
    ref_universe = mda.Universe(reference)
    analysis = rms.RMSD(universe, ref_universe, select=selection)
    analysis.run()
    output.write_text(json.dumps(analysis.rmsd.tolist(), indent=2), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("topology", type=Path)
    parser.add_argument("trajectory", type=Path)
    parser.add_argument("--occupancy", help="Atom selection for occupancy map")
    parser.add_argument("--contacts", nargs=3, metavar=("SEL1", "SEL2", "RADIUS"))
    parser.add_argument("--rmsd", nargs=2, metavar=("SELECTION", "REFERENCE"))
    parser.add_argument("--outdir", type=Path, default=Path("../analysis"))
    args = parser.parse_args()

    universe = mda.Universe(str(args.topology), str(args.trajectory))
    args.outdir.mkdir(parents=True, exist_ok=True)

    if args.occupancy:
        occupancy(universe, args.occupancy, args.outdir / "occupancy.dx")
    if args.contacts:
        sel1, sel2, radius = args.contacts
        contact_frequency(universe, sel1, sel2, float(radius), args.outdir / "contacts.json")
    if args.rmsd:
        selection, reference = args.rmsd
        compute_rmsd(universe, selection, reference, args.outdir / "rmsd.json")


if __name__ == "__main__":
    main()
