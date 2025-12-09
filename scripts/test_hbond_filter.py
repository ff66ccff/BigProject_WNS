#!/usr/bin/env python3
"""Test script for the stable hydrogen bond filter."""

from __future__ import annotations

import os
import sys
from pathlib import Path

# Add scripts directory to path
SCRIPTS_DIR = Path(__file__).resolve().parent
sys.path.append(str(SCRIPTS_DIR))

try:
    from filter_stable_hbonds import filter_and_save_hbond_complex
except ImportError as exc:
    print(f"❌ 错误: 无法导入filter_stable_hbonds模块 - {exc}")
    sys.exit(1)


def test_filter():
    """Test the hydrogen bond filter with sample data."""
    print("🧪 测试稳态氢键筛选器...")
    
    # Check if we have test data
    test_tpr = "gmx/final.tpr"
    test_traj = "gmx/final.xtc"
    test_output = "test_results/stable_hbond_complex.pdb"
    
    if not os.path.exists(test_tpr):
        print(f"⚠️ 警告: 测试文件不存在: {test_tpr}")
        print("请先运行完整的WnS流程生成测试数据")
        return False
    
    # Run the filter
    success = filter_and_save_hbond_complex(
        tpr_file=test_tpr,
        traj_file=test_traj if os.path.exists(test_traj) else None,
        output_pdb=test_output
    )
    
    if success:
        print("✅ 测试成功!")
        print(f"📄 结果文件: {test_output}")
        
        # Verify the output file
        if os.path.exists(test_output):
            print(f"✅ 输出文件已生成: {test_output}")
            
            # Quick check of the file content
            with open(test_output, 'r') as f:
                lines = f.readlines()
                atom_lines = [line for line in lines if line.startswith('ATOM')]
                ligand_lines = [line for line in atom_lines if line[17:20].strip() == 'LIG']
                
                print(f"📊 输出统计:")
                print(f"  - 总原子数: {len(atom_lines)}")
                print(f"  - 配体原子数: {len(ligand_lines)}")
                
                if ligand_lines:
                    ligand_resids = set()
                    for line in ligand_lines:
                        resid = int(line[22:26])
                        ligand_resids.add(resid)
                    
                    print(f"  - 配体残基数: {len(ligand_resids)}")
                    print(f"  - 配体残基ID: {sorted(ligand_resids)}")
        else:
            print("❌ 错误: 输出文件未生成")
            return False
    else:
        print("❌ 测试失败")
        return False
    
    return True


def main():
    """Main function."""
    print("🔬 稳态氢键筛选器测试")
    print("=" * 50)
    
    success = test_filter()
    
    if success:
        print("\n🎉 所有测试通过!")
        print("💡 提示: 稳态氢键筛选器工作正常")
    else:
        print("\n❌ 测试失败")
        print("💡 提示: 请检查输入数据或配置")
        sys.exit(1)


if __name__ == "__main__":
    main()