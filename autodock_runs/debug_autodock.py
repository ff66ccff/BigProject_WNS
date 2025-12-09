#!/usr/bin/env python3
"""调试AutoDock的map期望"""

import sys
import os

def analyze_dpf_file(dpf_path):
    """分析DPF文件的内容"""
    print(f"Analyzing DPF file: {dpf_path}")
    
    with open(dpf_path, 'r') as f:
        lines = f.readlines()
    
    ligand_types = None
    map_count = 0
    
    for line in lines:
        line = line.strip()
        if line.startswith('ligand_types'):
            ligand_types = line.split()[1:]
            print(f"Found ligand_types: {ligand_types}")
            print(f"Number of ligand types: {len(ligand_types)}")
        elif line.startswith('map'):
            map_count += 1
            print(f"Found map: {line}")
    
    print(f"\nSummary:")
    print(f"  ligand_types: {ligand_types}")
    print(f"  number of ligand types: {len(ligand_types) if ligand_types else 0}")
    print(f"  map count: {map_count}")
    
    if ligand_types and map_count == len(ligand_types):
        print("  ✓ DPF file is consistent")
        return True
    else:
        print("  ✗ DPF file has inconsistency")
        return False

if __name__ == "__main__":
    analyze_dpf_file("wrapper_101.dpf")