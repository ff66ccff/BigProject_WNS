#!/usr/bin/env python3
"""
修复配体拓扑文件中的二面角参数问题
"""

import sys
from pathlib import Path

def fix_topology_file(input_file, output_file):
    """修复拓扑文件中的二面角参数"""
    with open(input_file, 'r') as f_in, open(output_file, 'w') as f_out:
        lines = f_in.readlines()
        
        in_dihedrals = False
        dihedral_count = 0
        
        for line in lines:
            if line.strip().startswith('[ dihedrals ]'):
                in_dihedrals = True
                f_out.write(line)
                f_out.write("; ai   aj   ak   al   funct   phi   force.c.\n")
                continue
            
            if in_dihedrals and line.strip() and not line.startswith(';'):
                # 解析二面角行
                parts = line.split()
                if len(parts) >= 5:
                    ai, aj, ak, al, funct = parts[:5]
                    
                    # 检查是否有重复的原子索引
                    atoms = [int(ai), int(aj), int(ak), int(al)]
                    if len(set(atoms)) < 4:
                        continue  # 跳过重复原子的二面角
                    
                    # 修复参数格式 - GROMACS需要6个参数或3个参数
                    # 我们使用简化的3参数格式：ai aj ak al funct
                    dihedral_count += 1
                    f_out.write(f"{ai:>5} {aj:>5} {ak:>5} {al:>5}     1\n")
                continue
            
            f_out.write(line)
    
    print(f"修复完成，生成了 {dihedral_count} 个二面角")

def main():
    if len(sys.argv) != 3:
        print("用法: python fix_ligand_topology.py input.itp output.itp")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    print(f"修复拓扑文件: {input_file} -> {output_file}")
    fix_topology_file(input_file, output_file)
    print("完成！")

if __name__ == "__main__":
    main()