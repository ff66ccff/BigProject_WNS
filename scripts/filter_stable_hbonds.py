#!/usr/bin/env python3
"""Filter stable hydrogen bond complexes after Shaker simulation."""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

try:
    import MDAnalysis as mda
    from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis
except ImportError as exc:
    raise SystemExit("MDAnalysis is required. Install it via 'pip install MDAnalysis'.") from exc


def filter_and_save_hbond_complex(
    tpr_file: str = "shaker_results/final.tpr", 
    traj_file: str = "shaker_results/final.xtc", 
    output_pdb: str = "final_results/stable_hbond_complex.pdb"
) -> bool:
    """Filter and save complexes with stable hydrogen bonds.
    
    Args:
        tpr_file: GROMACS topology file
        traj_file: GROMACS trajectory file (optional)
        output_pdb: Output PDB file path
        
    Returns:
        True if successful, False if no ligands passed the filter
    """
    print(f"🔍 开始筛选稳态氢键配体...")
    print(f"拓扑文件: {tpr_file}")
    if traj_file and os.path.exists(traj_file):
        print(f"轨迹文件: {traj_file}")
    print(f"输出文件: {output_pdb}")
    
    # Check if input files exist
    if not os.path.exists(tpr_file):
        print(f"❌ 错误: 拓扑文件不存在: {tpr_file}")
        return False
    
    try:
        if traj_file and os.path.exists(traj_file):
            u = mda.Universe(tpr_file, traj_file)
            print(f"✅ 成功加载拓扑和轨迹文件")
        else:
            u = mda.Universe(tpr_file)
            print(f"✅ 成功加载拓扑文件")
    except Exception as e:
        print(f"❌ 错误: 无法加载文件 - {e}")
        return False

    # 1. 切换到最后一帧 (代表经过了高温震荡后的最终状态)
    if len(u.trajectory) > 1:
        u.trajectory[-1]
        print(f"✅ 切换到最后一帧 (帧 {u.trajectory.frame})")
    else:
        print(f"✅ 使用单帧结构")
    
    # 2. 定义氢键分析器
    # distance=3.5A, angle=120度 是通用标准
    hbonds = HydrogenBondAnalysis(
        universe=u,
        donors_selection="protein or resname LIG",
        acceptors_selection="protein or resname LIG",
        distance=3.5,
        angle=120.0
    )
    
    # 运行分析 (只分析当前这一帧)
    print("🔬 运行氢键分析...")
    hbonds.run(start=-1, stop=None)
    
    # 3. 筛选逻辑：找出哪些配体形成了氢键
    ligand_atoms = u.select_atoms("resname LIG")
    ligand_indices = set(ligand_atoms.indices)
    
    valid_ligand_resids = set()
    
    # 遍历所有检测到的氢键
    if hbonds.results.hbonds is not None and len(hbonds.results.hbonds) > 0:
        print(f"📊 检测到 {len(hbonds.results.hbonds)} 个氢键")
        for bond in hbonds.results.hbonds:
            # bond 格式: [frame, donor_idx, acceptor_idx, dist, angle]
            d_idx = int(bond[1])
            a_idx = int(bond[2])
            
            # 检查：如果 Donor 或 Acceptor 是配体的一部分，记录该配体的 Residue ID
            if d_idx in ligand_indices:
                atom = u.atoms[d_idx]
                valid_ligand_resids.add(atom.resid)
            elif a_idx in ligand_indices:
                atom = u.atoms[a_idx]
                valid_ligand_resids.add(atom.resid)
    else:
        print("⚠️ 警告: 没有检测到氢键")
    
    total_ligands = len(set(ligand_atoms.resids))
    valid_ligands = len(valid_ligand_resids)
    
    print(f"📊 检测结果: 总配体数 {total_ligands} -> 含有氢键的配体数 {valid_ligands}")
    
    # 4. 保存结果
    if not valid_ligand_resids:
        print("⚠️ 警告: 没有配体通过氢键筛选。")
        return False
    
    # 构建选择语句：只选蛋白 + 通过筛选的配体
    # 格式: "protein or (resname LIG and resid 1 5 9)"
    resid_str = " ".join(str(r) for r in sorted(valid_ligand_resids))
    selection_cmd = f"protein or (resname LIG and resid {resid_str})"
    
    print(f"🎯 选择语句: {selection_cmd}")
    
    final_system = u.select_atoms(selection_cmd)
    
    # 确保输出目录存在
    output_dir = os.path.dirname(output_pdb)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
        print(f"📁 创建输出目录: {output_dir}")
    
    final_system.write(output_pdb)
    print(f"✅ 成功保存文件: {output_pdb}")
    print(f"📈 保存了 {len(final_system)} 个原子")
    
    # 打印通过的配体详细信息
    print("\n📋 通过筛选的配体详情:")
    for resid in sorted(valid_ligand_resids):
        lig_atoms = u.select_atoms(f"resname LIG and resid {resid}")
        print(f"  - 配体 {resid}: {len(lig_atoms)} 个原子")
    
    return True


def main() -> None:
    """Main function with argument parsing."""
    parser = argparse.ArgumentParser(
        description="Filter complexes with stable hydrogen bonds after Shaker simulation"
    )
    parser.add_argument(
        "--tpr", 
        default="shaker_results/final.tpr",
        help="GROMACS topology file (default: shaker_results/final.tpr)"
    )
    parser.add_argument(
        "--traj", 
        default="shaker_results/final.xtc",
        help="GROMACS trajectory file (default: shaker_results/final.xtc)"
    )
    parser.add_argument(
        "--output", 
        default="final_results/stable_hbond_complex.pdb",
        help="Output PDB file (default: final_results/stable_hbond_complex.pdb)"
    )
    parser.add_argument(
        "--distance", 
        type=float, 
        default=3.5,
        help="Hydrogen bond distance cutoff in Angstroms (default: 3.5)"
    )
    parser.add_argument(
        "--angle", 
        type=float, 
        default=120.0,
        help="Hydrogen bond angle cutoff in degrees (default: 120.0)"
    )
    
    args = parser.parse_args()
    
    # Convert to Path objects for consistency
    tpr_file = Path(args.tpr)
    traj_file = Path(args.traj)
    output_pdb = Path(args.output)
    
    # Check if topology file exists
    if not tpr_file.exists():
        print(f"❌ 错误: 拓扑文件不存在: {tpr_file}")
        sys.exit(1)
    
    # Check if trajectory file exists (optional)
    if traj_file.exists():
        print(f"📁 使用轨迹文件: {traj_file}")
        traj_file_str = str(traj_file)
    else:
        print(f"⚠️ 轨迹文件不存在，仅使用拓扑文件: {traj_file}")
        traj_file_str = None
    
    success = filter_and_save_hbond_complex(
        tpr_file=str(tpr_file),
        traj_file=traj_file_str,
        output_pdb=str(output_pdb)
    )
    
    if not success:
        print("❌ 筛选失败: 没有配体通过氢键筛选")
        sys.exit(1)
    
    print("\n🎉 稳态氢键筛选完成!")
    print(f"📄 结果文件: {output_pdb}")
    print("💡 提示: 该文件中的每个配体都经过了高温震荡并保持至少一个氢键")


if __name__ == "__main__":
    main()