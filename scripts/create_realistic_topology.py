#!/usr/bin/env python3
"""基于实际配体坐标生成真实的拓扑文件"""

import sys
import math
from pathlib import Path

def read_gro_file(gro_file):
    """读取GRO文件中的配体原子信息"""
    atoms = []
    with open(gro_file, 'r') as f:
        lines = f.readlines()
    
    # 找到配体原子（通常在蛋白质后面）
    for line in lines:
        if line.strip() and not line.startswith(';'):
            parts = line.split()
            if len(parts) >= 5:
                resname = parts[1]
                if resname == 'UNL':  # 配体残基名
                    atom_id = int(parts[2])
                    atom_name = parts[3]
                    x = float(parts[4])
                    y = float(parts[5])
                    z = float(parts[6])
                    atoms.append({
                        'id': atom_id,
                        'name': atom_name,
                        'x': x, 'y': y, 'z': z,
                        'type': atom_name[0]  # 使用原子名的第一个字符作为类型
                    })
    return atoms

def determine_bonds(atoms, cutoff_factor=1.2):
    """基于原子间距离确定键连接"""
    bonds = []
    
    # 原子半径（用于确定键长阈值）
    atomic_radii = {
        'H': 0.31, 'C': 0.76, 'N': 0.71, 'O': 0.66,
        'S': 1.05, 'P': 1.07
    }
    
    # 典型键长
    bond_lengths = {
        'H-H': 0.74, 'C-C': 1.54, 'C-H': 1.09, 'C-N': 1.47,
        'C-O': 1.43, 'N-H': 1.01, 'O-H': 0.96, 'N-O': 1.40,
        'C=S': 1.56, 'C=P': 1.84, 'S-H': 1.34, 'P-H': 1.42
    }
    
    for i in range(len(atoms)):
        for j in range(i+1, len(atoms)):
            atom1 = atoms[i]
            atom2 = atoms[j]
            
            # 计算距离
            dx = atom1['x'] - atom2['x']
            dy = atom1['y'] - atom2['y']
            dz = atom1['z'] - atom2['z']
            distance = math.sqrt(dx*dx + dy*dy + dz*dz)
            
            # 确定键类型
            type1 = atom1['type']
            type2 = atom2['type']
            bond_type = f"{type1}-{type2}" if type1 <= type2 else f"{type2}-{type1}"
            
            # 获取典型键长
            if bond_type in bond_lengths:
                typical_length = bond_lengths[bond_type]
            else:
                # 使用原子半径估算
                r1 = atomic_radii.get(type1, 0.77)
                r2 = atomic_radii.get(type2, 0.77)
                typical_length = (r1 + r2) * 1.1
            
            # 如果距离在合理范围内，认为是键
            if distance < typical_length * cutoff_factor:
                bonds.append({
                    'ai': atom1['id'],
                    'aj': atom2['id'],
                    'type': bond_type,
                    'length': typical_length,
                    'distance': distance
                })
    
    return bonds

def get_bond_parameters(bond_type):
    """获取键的力常数参数"""
    force_constants = {
        'H-H': 500, 'C-C': 350, 'C-H': 450, 'C-N': 380,
        'C-O': 420, 'N-H': 500, 'O-H': 550, 'N-O': 400,
        'C=S': 300, 'C=P': 280, 'S-H': 350, 'P-H': 400
    }
    return force_constants.get(bond_type, 400)

def create_topology(atoms, bonds):
    """创建拓扑文件"""
    topology = """; 基于实际结构的配体拓扑文件
; 使用真实的键连接和参数

[ moleculetype ]
; molname   nrexcl
UNL         3

[ atoms ]
; nr  type  resnr  residu atom  cgnr  charge   mass
"""
    
    # 添加原子定义
    atomic_masses = {'H': 1.008, 'C': 12.011, 'N': 14.007, 'O': 15.999, 'S': 32.066, 'P': 30.974}
    
    for atom in atoms:
        atom_type = atom['type']
        mass = atomic_masses.get(atom_type, 12.011)
        topology += f"  {atom['id']:4d}  {atom_type:5s}  1 UNL    {atom['name']:4s}  {atom['id']:4d}     0.000  {mass:6.3f}\n"
    
    topology += """
[ bonds ]
; ai   aj   funct   length   force.c.
"""
    
    # 添加键定义
    for bond in bonds:
        force_const = get_bond_parameters(bond['type'])
        topology += f"  {bond['ai']:4d}  {bond['aj']:4d}     1    {bond['length']:6.3f}    {force_const:6.0f}\n"
    
    return topology

def main():
    if len(sys.argv) != 3:
        print("用法: python create_realistic_topology.py <gro_file> <output_itp>")
        sys.exit(1)
    
    gro_file = sys.argv[1]
    output_file = sys.argv[2]
    
    print(f"读取GRO文件: {gro_file}")
    atoms = read_gro_file(gro_file)
    print(f"找到 {len(atoms)} 个配体原子")
    
    print("确定键连接...")
    bonds = determine_bonds(atoms)
    print(f"找到 {len(bonds)} 个键")
    
    print("生成拓扑文件...")
    topology = create_topology(atoms, bonds)
    
    with open(output_file, 'w') as f:
        f.write(topology)
    
    print(f"拓扑文件已保存到: {output_file}")
    
    # 统计键类型
    bond_types = {}
    for bond in bonds:
        bond_type = bond['type']
        bond_types[bond_type] = bond_types.get(bond_type, 0) + 1
    
    print("键类型统计:")
    for bond_type, count in sorted(bond_types.items()):
        print(f"  {bond_type}: {count}")

if __name__ == "__main__":
    main()