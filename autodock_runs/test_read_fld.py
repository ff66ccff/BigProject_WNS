#!/usr/bin/env python3
"""测试AutoDock是否能正确读取FLD文件"""

import sys
import os

def check_fld_file(fld_path):
    """检查FLD文件的内容"""
    print(f"Checking FLD file: {fld_path}")
    
    with open(fld_path, 'r') as f:
        lines = f.readlines()
    
    veclen = None
    variables = []
    labels = []
    
    for line in lines:
        line = line.strip()
        if line.startswith('veclen='):
            veclen = int(line.split('=')[1].split()[0])
            print(f"Found veclen={veclen}")
        elif line.startswith('variable '):
            parts = line.split()
            if len(parts) >= 4:
                var_num = int(parts[1])
                filename = parts[3]
                variables.append((var_num, filename))
                print(f"Found variable {var_num}: {filename}")
        elif line.startswith('label='):
            label = line.split('=')[1].strip()
            labels.append(label)
            print(f"Found label: {label}")
    
    print(f"\nSummary:")
    print(f"  veclen: {veclen}")
    print(f"  variables found: {len(variables)}")
    print(f"  labels found: {len(labels)}")
    
    if veclen and len(variables) == veclen:
        print("  ✓ FLD file is consistent")
        return True
    else:
        print("  ✗ FLD file has inconsistency")
        return False

if __name__ == "__main__":
    check_fld_file("protein.maps.fld")