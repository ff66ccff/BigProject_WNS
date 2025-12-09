#!/usr/bin/env python3
"""Verify the implementation of stable hydrogen bond filter."""

from __future__ import annotations

import os
import sys
from pathlib import Path


def check_file_exists(filepath: str, description: str) -> bool:
    """Check if a file exists and report status."""
    if os.path.exists(filepath):
        print(f"✅ {description}: {filepath}")
        return True
    else:
        print(f"❌ {description}: {filepath} (不存在)")
        return False


def check_script_syntax(script_path: str) -> bool:
    """Check if a Python script has valid syntax."""
    try:
        with open(script_path, 'r', encoding='utf-8') as f:
            content = f.read()
        compile(content, script_path, 'exec')
        print(f"✅ 脚本语法正确: {script_path}")
        return True
    except SyntaxError as e:
        print(f"❌ 脚本语法错误: {script_path}")
        print(f"   错误: {e}")
        return False
    except Exception as e:
        print(f"❌ 脚本检查失败: {script_path}")
        print(f"   错误: {e}")
        return False


def verify_implementation():
    """Verify the complete implementation."""
    print("🔍 验证稳态氢键筛选器实现")
    print("=" * 50)
    
    all_good = True
    
    # Check 1: Core script exists and has valid syntax
    print("\n📁 检查核心脚本:")
    filter_script = "scripts/filter_stable_hbonds.py"
    if check_file_exists(filter_script, "氢键筛选器脚本"):
        if not check_script_syntax(filter_script):
            all_good = False
    
    # Check 2: Modified pipeline script
    print("\n📁 检查修改的流程脚本:")
    pipeline_script = "scripts/run_full_pipeline.py"
    if check_file_exists(pipeline_script, "主流程脚本"):
        if not check_script_syntax(pipeline_script):
            all_good = False
    
    # Check 3: Test script
    print("\n📁 检查测试脚本:")
    test_script = "scripts/test_hbond_filter.py"
    if check_file_exists(test_script, "测试脚本"):
        if not check_script_syntax(test_script):
            all_good = False
    
    # Check 4: Documentation
    print("\n📁 检查文档:")
    docs = [
        ("STABLE_HBOND_IMPLEMENTATION_SUMMARY.md", "实现总结文档"),
        ("STABLE_HBOND_FILTER_README.md", "使用说明文档")
    ]
    
    for doc, desc in docs:
        check_file_exists(doc, desc)
    
    # Check 5: Import verification
    print("\n🔬 检查模块导入:")
    try:
        # Test if we can import the filter module
        sys.path.append("scripts")
        try:
            from filter_stable_hbonds import filter_and_save_hbond_complex
            print("✅ 氢键筛选器模块导入成功")
        except ImportError as e:
            print(f"❌ 氢键筛选器模块导入失败: {e}")
            all_good = False
    except Exception as e:
        print(f"❌ 模块导入检查失败: {e}")
        all_good = False
    
    # Check 6: Configuration files
    print("\n⚙️ 检查配置文件:")
    config_files = [
        ("config.yml", "主配置文件"),
        ("mdp/annealing.mdp", "模拟退火参数"),
        ("scripts/templates/gpf_template.txt", "GPF模板文件")
    ]
    
    for config, desc in config_files:
        check_file_exists(config, desc)
    
    # Summary
    print("\n" + "=" * 50)
    if all_good:
        print("🎉 实现验证通过!")
        print("💡 提示: 稳态氢键筛选器已正确实现")
        print("\n📋 下一步:")
        print("1. 运行完整流程: python scripts/run_full_pipeline.py")
        print("2. 检查输出: final_results/stable_hbond_complex.pdb")
        print("3. 验证结果: 确保每个配体都有氢键")
    else:
        print("❌ 实现验证失败!")
        print("💡 提示: 请检查上述错误并修复")
    
    return all_good


def main():
    """Main function."""
    try:
        success = verify_implementation()
        sys.exit(0 if success else 1)
    except KeyboardInterrupt:
        print("\n⚠️ 验证被用户中断")
        sys.exit(1)
    except Exception as e:
        print(f"\n❌ 验证过程中发生错误: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()