#!/usr/bin/env python3
"""
穷人的拓扑生成器 - 为配体创建基本的键合参数
基于距离和化学常识生成键、角、二面角参数
"""

import sys
import math
from pathlib import Path

def read_ligand_coordinates(gro_file):
    """读取GRO文件中的配体坐标"""
    coords = {}
    atoms = []
    
    with open(gro_file, 'r') as f:
        lines = f.readlines()
    
    # 找到配体原子（假设配体残基名为UNL）
    in_ligand = False
    atom_id = 1
    
    for line in lines[2:]:  # 跳过标题行
        if len(line) < 40:
            continue
            
        resname = line[5:10].strip()
        if resname == 'UNL':
            in_ligand = True
            atom_name = line[10:15].strip()
            x = float(line[20:28])
            y = float(line[28:36])
            z = float(line[36:44])
            
            atoms.append({
                'id': atom_id,
                'name': atom_name,
                'resname': resname,
                'x': x, 'y': y, 'z': z
            })
            coords[atom_id] = (x, y, z)
            atom_id += 1
        elif in_ligand and resname != 'UNL':
            break
    
    return atoms, coords

def calculate_distance(coord1, coord2):
    """计算两点间距离"""
    dx = coord1[0] - coord2[0]
    dy = coord1[1] - coord2[1]
    dz = coord1[2] - coord2[2]
    return math.sqrt(dx*dx + dy*dy + dz*dz)

def determine_bond_type(atom1_name, atom2_name, distance):
    """根据原子类型和距离确定键类型"""
    # 基本键长阈值（Å）
    single_bond_max = {
        ('C', 'C'): 1.6,
        ('C', 'H'): 1.2,
        ('C', 'N'): 1.5,
        ('C', 'O'): 1.5,
        ('C', 'S'): 1.9,
        ('N', 'H'): 1.1,
        ('N', 'N'): 1.5,
        ('N', 'O'): 1.5,
        ('O', 'H'): 1.1,
        ('O', 'O'): 1.5,
        ('S', 'H'): 1.4,
        ('S', 'S'): 2.2
    }
    
    # 获取原子类型（第一个字符）
    type1 = atom1_name[0] if atom1_name else 'C'
    type2 = atom2_name[0] if atom2_name else 'C'
    
    # 检查键类型
    key = (type1, type2)
    if key not in single_bond_max:
        key = (type2, type1)
    
    if key in single_bond_max:
        max_dist = single_bond_max[key]
        if distance <= max_dist:
            return 'single'
        elif distance <= max_dist * 1.3:
            return 'double'
        else:
            return None
    else:
        # 默认判断
        if distance <= 1.6:
            return 'single'
        elif distance <= 2.0:
            return 'double'
        else:
            return None

def generate_bonds(atoms, coords):
    """生成键列表"""
    bonds = []
    bond_id = 1
    
    for i, atom1 in enumerate(atoms):
        for atom2 in atoms[i+1:]:
            dist = calculate_distance(coords[atom1['id']], coords[atom2['id']])
            bond_type = determine_bond_type(atom1['name'], atom2['name'], dist)
            
            if bond_type:
                bonds.append({
                    'id': bond_id,
                    'ai': atom1['id'],
                    'aj': atom2['id'],
                    'type': bond_type
                })
                bond_id += 1
    
    return bonds

def generate_angles(bonds):
    """生成角列表"""
    angles = []
    angle_id = 1
    
    # 构建键邻接表
    bond_map = {}
    for bond in bonds:
        if bond['ai'] not in bond_map:
            bond_map[bond['ai']] = []
        if bond['aj'] not in bond_map:
            bond_map[bond['aj']] = []
        bond_map[bond['ai']].append(bond['aj'])
        bond_map[bond['aj']].append(bond['ai'])
    
    # 生成角
    for center_atom, neighbors in bond_map.items():
        if len(neighbors) >= 2:
            for i in range(len(neighbors)):
                for j in range(i+1, len(neighbors)):
                    angles.append({
                        'id': angle_id,
                        'ai': neighbors[i],
                        'aj': center_atom,
                        'ak': neighbors[j],
                        'type': 'normal'
                    })
                    angle_id += 1
    
    return angles

def generate_dihedrals(bonds):
    """生成二面角列表"""
    dihedrals = []
    dihedral_id = 1
    
    # 构建键邻接表
    bond_map = {}
    for bond in bonds:
        if bond['ai'] not in bond_map:
            bond_map[bond['ai']] = []
        if bond['aj'] not in bond_map:
            bond_map[bond['aj']] = []
        bond_map[bond['ai']].append(bond['aj'])
        bond_map[bond['aj']].append(bond['ai'])
    
    # 生成二面角
    for bond in bonds:
        ai, aj = bond['ai'], bond['aj']
        
        # 找到ai的其他邻居
        if ai in bond_map and len(bond_map[ai]) > 1:
            for ak in bond_map[ai]:
                if ak != aj:
                    # 找到aj的其他邻居
                    if aj in bond_map and len(bond_map[aj]) > 1:
                        for al in bond_map[aj]:
                            if al != ai:
                                dihedrals.append({
                                    'id': dihedral_id,
                                    'ai': ak,
                                    'aj': ai,
                                    'ak': aj,
                                    'al': al,
                                    'type': 'normal'
                                })
                                dihedral_id += 1
    
    return dihedrals

def write_topology_file(atoms, bonds, angles, dihedrals, output_file):
    """写入拓扑文件"""
    with open(output_file, 'w') as f:
        f.write("; 穷人的拓扑生成器 - 自动生成的配体拓扑\n")
        f.write("; 警告：这是基于距离的简单估算，可能不准确\n")
        f.write("; 建议使用量子化学计算或实验数据验证\n\n")
        
        f.write("[ moleculetype ]\n")
        f.write("; molname   nrexcl\n")
        f.write("UNL         3\n\n")
        
        f.write("[ atoms ]\n")
        f.write("; nr  type  resnr  residu atom  cgnr  charge   mass\n")
        for atom in atoms:
            # 简化的原子类型分配
            atom_type = atom['name'][0]  # 使用第一个字符作为类型
            if atom_type == 'C':
                gromacs_type = 'C'
                mass = 12.011
                charge = 0.0
            elif atom_type == 'H':
                gromacs_type = 'H'
                mass = 1.008
                charge = 0.0
            elif atom_type == 'N':
                gromacs_type = 'N'
                mass = 14.007
                charge = 0.0
            elif atom_type == 'O':
                gromacs_type = 'O'
                mass = 15.999
                charge = 0.0
            elif atom_type == 'S':
                gromacs_type = 'S'
                mass = 32.066
                charge = 0.0
            else:
                gromacs_type = 'C'  # 默认
                mass = 12.011
                charge = 0.0
            
            f.write(f"{atom['id']:>5} {gromacs_type:>5}     1 UNL     {atom['name']:>3} {atom['id']:>5}   {charge:>7.3f}  {mass:>6.3f}\n")
        
        f.write("\n[ bonds ]\n")
        f.write("; ai   aj   funct   length   force.c.\n")
        for bond in bonds:
            # 简化的键参数
            if bond['type'] == 'single':
                length = 1.5  # nm
                force = 400   # kJ/mol/nm^2
            elif bond['type'] == 'double':
                length = 1.3  # nm
                force = 600   # kJ/mol/nm^2
            else:
                length = 1.4  # nm
                force = 500   # kJ/mol/nm^2
            
            f.write(f"{bond['ai']:>5} {bond['aj']:>5}     1   {length:>6.3f}   {force:>6.0f}\n")
        
        f.write("\n[ angles ]\n")
        f.write("; ai   aj   ak   funct   angle   force.c.\n")
        for angle in angles:
            # 简化的角参数
            angle_deg = 109.5  # degrees
            force = 300       # kJ/mol/rad^2
            f.write(f"{angle['ai']:>5} {angle['aj']:>5} {angle['ak']:>5}     1   {angle_deg:>6.1f}   {force:>6.0f}\n")
        
        f.write("\n[ dihedrals ]\n")
        f.write("; ai   aj   ak   al   funct   phi   force.c.\n")
        for dihedral in dihedrals:
            # 简化的二面角参数
            phi = 180.0    # degrees
            force = 5       # kJ/mol
            f.write(f"{dihedral['ai']:>5} {dihedral['aj']:>5} {dihedral['ak']:>5} {dihedral['al']:>5}     1   {phi:>6.1f}   {force:>6.0f}\n")

def main():
    if len(sys.argv) != 3:
        print("用法: python poor_mans_topology.py input.gro output.itp")
        sys.exit(1)
    
    input_gro = sys.argv[1]
    output_itp = sys.argv[2]
    
    print(f"读取配体坐标: {input_gro}")
    atoms, coords = read_ligand_coordinates(input_gro)
    print(f"找到 {len(atoms)} 个配体原子")
    
    print("生成键...")
    bonds = generate_bonds(atoms, coords)
    print(f"生成 {len(bonds)} 个键")
    
    print("生成角...")
    angles = generate_angles(bonds)
    print(f"生成 {len(angles)} 个角")
    
    print("生成二面角...")
    dihedrals = generate_dihedrals(bonds)
    print(f"生成 {len(dihedrals)} 个二面角")
    
    print(f"写入拓扑文件: {output_itp}")
    write_topology_file(atoms, bonds, angles, dihedrals, output_itp)
    
    print("完成！")

if __name__ == "__main__":
    main()