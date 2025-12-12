#!/usr/bin/env python3
"""
简化版配体拓扑生成器 - 只处理配体分子
"""

import sys
import math
from pathlib import Path

def extract_ligand_gro(input_gro, output_ligand_gro):
    """从完整系统中提取配体分子"""
    with open(input_gro, 'r') as f_in, open(output_ligand_gro, 'w') as f_out:
        lines = f_in.readlines()
        
        # 写入标题
        f_out.write("Ligand extracted from system\n")
        f_out.write("568\n")  # 配体原子数
        
        # 提取配体原子
        atom_id = 1
        for line in lines[2:]:  # 跳过标题行
            if len(line) < 40:
                continue
                
            resname = line[5:10].strip()
            if resname == 'UNL':
                # 重写原子ID
                atom_name = line[10:15].strip()
                x = float(line[20:28])
                y = float(line[28:36])
                z = float(line[36:44])
                
                f_out.write(f"    1UNL{atom_name:>5} {atom_id:>5} {x:>8.3f}{y:>8.3f}{z:>8.3f}\n")
                atom_id += 1
            elif atom_id > 1:  # 已经开始处理配体，现在结束了
                break
        
        # 写入盒子尺寸
        f_out.write(lines[-1])

def read_ligand_coordinates(gro_file):
    """读取配体坐标"""
    coords = {}
    atoms = []
    
    with open(gro_file, 'r') as f:
        lines = f.readlines()
    
    atom_id = 1
    for line in lines[2:-1]:  # 跳过标题和盒子尺寸
        if len(line) < 40:
            continue
            
        atom_name = line[10:15].strip()
        x = float(line[20:28])
        y = float(line[28:36])
        z = float(line[36:44])
        
        atoms.append({
            'id': atom_id,
            'name': atom_name,
            'x': x, 'y': y, 'z': z
        })
        coords[atom_id] = (x, y, z)
        atom_id += 1
    
    return atoms, coords

def calculate_distance(coord1, coord2):
    """计算两点间距离"""
    dx = coord1[0] - coord2[0]
    dy = coord1[1] - coord2[1]
    dz = coord1[2] - coord2[2]
    return math.sqrt(dx*dx + dy*dy + dz*dz)

def generate_bonds(atoms, coords):
    """生成键列表 - 只考虑最近邻"""
    bonds = []
    bond_id = 1
    
    for i, atom1 in enumerate(atoms):
        distances = []
        for j, atom2 in enumerate(atoms):
            if i != j:
                dist = calculate_distance(coords[atom1['id']], coords[atom2['id']])
                distances.append((dist, atom2['id']))
        
        # 排序并取最近的几个原子作为键合对象
        distances.sort()
        
        # 根据原子类型决定键的数量
        atom1_type = atom1['name'][0]
        max_bonds = {'C': 4, 'N': 3, 'O': 2, 'S': 2, 'H': 1}.get(atom1_type, 4)
        
        bond_count = 0
        for dist, atom2_id in distances:
            if bond_count >= max_bonds:
                break
                
            # 距离阈值
            if dist <= 1.6:  # 单键阈值
                # 检查是否已经存在这个键
                exists = False
                for bond in bonds:
                    if (bond['ai'] == atom1['id'] and bond['aj'] == atom2_id) or \
                       (bond['ai'] == atom2_id and bond['aj'] == atom1['id']):
                        exists = True
                        break
                
                if not exists:
                    bonds.append({
                        'id': bond_id,
                        'ai': atom1['id'],
                        'aj': atom2_id,
                        'type': 'single'
                    })
                    bond_id += 1
                    bond_count += 1
    
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
                                break  # 只取一个，避免过多
                        break  # 只取一个，避免过多
    
    return dihedrals

def write_topology_file(atoms, bonds, angles, dihedrals, output_file):
    """写入拓扑文件"""
    with open(output_file, 'w') as f:
        f.write("; 简化版配体拓扑生成器\n")
        f.write("; 警告：这是基于距离的估算，可能不准确\n\n")
        
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
            length = 1.5  # nm
            force = 400   # kJ/mol/nm^2
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
        print("用法: python simple_ligand_topology.py input.gro output.itp")
        sys.exit(1)
    
    input_gro = sys.argv[1]
    output_itp = sys.argv[2]
    
    # 首先提取配体
    ligand_gro = "temp_ligand.gro"
    print(f"提取配体到: {ligand_gro}")
    extract_ligand_gro(input_gro, ligand_gro)
    
    print(f"读取配体坐标: {ligand_gro}")
    atoms, coords = read_ligand_coordinates(ligand_gro)
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
    
    # 清理临时文件
    Path(ligand_gro).unlink()
    
    print("完成！")

if __name__ == "__main__":
    main()