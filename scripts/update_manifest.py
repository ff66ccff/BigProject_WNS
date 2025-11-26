#!/usr/bin/env python3
"""Utility for updating the project manifest with reproducibility metadata."""

from __future__ import annotations

import argparse
import os
import sys
from typing import Any, Dict, List

try:
    import yaml  # type: ignore
except ImportError as exc:  # pragma: no cover
    raise SystemExit(
        "PyYAML is required to run update_manifest.py. Install it via 'pip install pyyaml'."
    ) from exc

DEFAULT_MANIFEST_PATH = os.path.join(
    os.path.dirname(os.path.dirname(__file__)), "manifest", "run-manifest.yml"
)


def load_manifest(path: str) -> Dict[str, Any]:
    if not os.path.exists(path):
        return {}
    with open(path, "r", encoding="utf-8") as handle:
        content = handle.read().strip()
        return yaml.safe_load(content) if content else {}


def save_manifest(path: str, manifest: Dict[str, Any]) -> None:
    with open(path, "w", encoding="utf-8") as handle:
        yaml.safe_dump(manifest, handle, sort_keys=False, allow_unicode=False)


def update_software(manifest: Dict[str, Any], items: List[str]) -> None:
    if "software" not in manifest:
        manifest["software"] = {}
    for item in items:
        try:
            key, value = item.split("=", 1)
        except ValueError:
            raise ValueError(f"Invalid software spec '{item}'. Use name=version format.")
        manifest["software"].setdefault(key, {})
        manifest["software"][key]["version"] = value


def update_kv_section(manifest: Dict[str, Any], section: str, items: List[str]) -> None:
    if section not in manifest or manifest[section] is None:
        manifest[section] = {}
    for item in items:
        try:
            key, value = item.split("=", 1)
        except ValueError:
            raise ValueError(f"Invalid {section} spec '{item}'. Use key=value format.")
        manifest[section][key] = value


def append_list(manifest: Dict[str, Any], section: str, values: List[Any]) -> None:
    manifest.setdefault(section, [])
    manifest[section].extend(values)


def update_nested_parameter(manifest: Dict[str, Any], entries: List[str]) -> None:
    manifest.setdefault("parameters", {})
    for entry in entries:
        try:
            qualified_key, value = entry.split("=", 1)
            section, key = qualified_key.split(".", 1)
        except ValueError:
            raise ValueError(
                f"Invalid parameter entry '{entry}'. Use section.key=value format."
            )
        manifest["parameters"].setdefault(section, {})
        manifest["parameters"][section][key] = value


def parse_args(argv: List[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--manifest",
        default=DEFAULT_MANIFEST_PATH,
        help="Path to the manifest YAML file (defaults to project/manifest/run-manifest.yml)",
    )
    parser.add_argument("--software", action="append", default=[], help="Add software name=version entry")
    parser.add_argument("--command", action="append", default=[], help="Append executed command string")
    parser.add_argument("--seed", action="append", default=[], help="Append integer seed value")
    parser.add_argument("--input", action="append", default=[], help="Add input key=path record")
    parser.add_argument("--output", action="append", default=[], help="Add output key=path record")
    parser.add_argument(
        "--parameter",
        action="append",
        default=[],
        help="Add parameter entry using section.key=value notation",
    )
    parser.add_argument(
        "--note",
        help="Replace the free-form notes field",
    )
    parser.add_argument(
        "--append-note",
        action="append",
        default=[],
        help="Append text to the free-form notes field",
    )
    return parser.parse_args(argv)


def main(argv: List[str]) -> None:
    args = parse_args(argv)
    manifest_path = os.path.abspath(args.manifest)
    manifest = load_manifest(manifest_path)

    try:
        if args.software:
            update_software(manifest, args.software)
        if args.command:
            append_list(manifest, "commands", args.command)
        if args.seed:
            append_list(manifest, "random_seeds", args.seed)
        if args.input:
            update_kv_section(manifest, "inputs", args.input)
        if args.output:
            update_kv_section(manifest, "outputs", args.output)
        if args.parameter:
            update_nested_parameter(manifest, args.parameter)
        if args.note is not None:
            manifest["notes"] = args.note
        if args.append_note:
            note_text = manifest.get("notes", "")
            note_lines = [note_text] if note_text else []
            note_lines.extend(args.append_note)
            manifest["notes"] = "\n".join(note_lines)
    except ValueError as exc:
        raise SystemExit(str(exc))

    save_manifest(manifest_path, manifest)


if __name__ == "__main__":
    main(sys.argv[1:])
