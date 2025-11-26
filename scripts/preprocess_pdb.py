#!/usr/bin/env python3
"""Minimal PDB cleaner: normalizes chain identifiers and removes duplicate atoms."""

from __future__ import annotations

import argparse
import itertools
import os
from typing import Iterable, List, Set, Tuple

CHAIN_IDS = list("ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789")


def _normalize_chain_ids(lines: Iterable[str]) -> List[str]:
    chain_iterator = (cid for cid in CHAIN_IDS)
    chain_map = {}
    normalized_lines: List[str] = []
    for line in lines:
        if line.startswith(("ATOM", "HETATM")) and len(line) >= 22:
            chain_id = line[21].strip() or ""
            residue_key = line[17:27]
            if chain_id == "":
                if residue_key not in chain_map:
                    try:
                        chain_map[residue_key] = next(chain_iterator)
                    except StopIteration:
                        raise RuntimeError(
                            "Ran out of chain identifiers while normalizing unnamed chains."
                        )
                chain_id = chain_map[residue_key]
            normalized_line = line[:21] + chain_id + line[22:]
            normalized_lines.append(normalized_line)
        else:
            normalized_lines.append(line)
    return normalized_lines


def _deduplicate_atoms(lines: Iterable[str]) -> List[str]:
    seen: Set[Tuple[str, str, str, str]] = set()
    deduped_lines: List[str] = []
    for line in lines:
        if line.startswith(("ATOM", "HETATM")) and len(line) >= 54:
            atom_name = line[12:16]
            res_id = line[17:27]
            chain_id = line[21]
            coords = line[30:54]
            key = (atom_name, res_id, chain_id, coords)
            if key in seen:
                continue
            seen.add(key)
            deduped_lines.append(line)
        else:
            deduped_lines.append(line)
    return deduped_lines


def clean_pdb(input_path: str, output_path: str) -> None:
    if not os.path.exists(input_path):
        raise FileNotFoundError(f"Input PDB '{input_path}' does not exist.")

    with open(input_path, "r", encoding="utf-8") as handle:
        original_lines = handle.readlines()

    normalized = _normalize_chain_ids(original_lines)
    cleaned = _deduplicate_atoms(normalized)

    with open(output_path, "w", encoding="utf-8") as handle:
        handle.writelines(cleaned)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", help="Path to the input PDB file")
    parser.add_argument(
        "-o",
        "--output",
        default="protein_clean.pdb",
        help="Destination path for the cleaned PDB (default: protein_clean.pdb)",
    )
    return parser


def main(argv: List[str] | None = None) -> None:
    parser = build_parser()
    args = parser.parse_args(argv)
    clean_pdb(args.input, args.output)


if __name__ == "__main__":
    main()
