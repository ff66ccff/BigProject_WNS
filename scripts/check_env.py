#!/usr/bin/env python3
"""Check environment for Wrap 'n' Shake MD simulation."""

import shutil
import sys
import subprocess
from pathlib import Path

def check_tool(name, package_suggestion=None):
    """Check if a tool is available in PATH."""
    path = shutil.which(name)
    status = "OK" if path else "MISSING"
    print(f"{name:15} : {status}")
    if not path and package_suggestion:
        print(f"  -> Install via: {package_suggestion}")
    return path is not None

def main():
    print("=== Checking MD Environment ===\n")
    
    # Check GROMACS
    has_gmx = check_tool("gmx", "sudo apt-get install gromacs")
    if has_gmx:
        # Check version
        try:
            output = subprocess.check_output(["gmx", "--version"], text=True)
            version = output.splitlines()[0] if output else "Unknown"
            print(f"  -> Version: {version}")
        except:
            pass

    # Check AmberTools components
    print("\n--- Ligand Parameterization Tools ---")
    has_acpype = check_tool("acpype", "pip install acpype (requires AmberTools or OpenBabel)")
    has_antechamber = check_tool("antechamber", "sudo apt-get install ambertools")
    has_parmchk2 = check_tool("parmchk2", "sudo apt-get install ambertools")
    
    print("\n--- Python Dependencies ---")
    try:
        import MDAnalysis
        print(f"{'MDAnalysis':15} : OK")
    except ImportError:
        print(f"{'MDAnalysis':15} : MISSING")
        print("  -> pip install MDAnalysis")

    print("\n=== Summary ===")
    if has_gmx and has_acpype:
        print("Environment looks good for basic MD!")
    else:
        print("Some dependencies are missing. Please install them before starting MD.")
        if not has_acpype:
             print("CRITICAL: acpype is required for ligand parameterization.")

if __name__ == "__main__":
    main()
