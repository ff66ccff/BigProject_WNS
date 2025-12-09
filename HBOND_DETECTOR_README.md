# 氢键探测器 (Hydrogen Bond Detector)

基于WnS (Wrap 'n' Shake) 算法的专门氢键探测工具，用于识别和可视化蛋白质-配体复合物中的稳定氢键。

## 🎯 功能概述

本工具通过三阶段流程识别稳定的氢键：

1. **Wrapper阶段**：在蛋白质表面生成多个配体构象
2. **Shaker阶段**：通过模拟退火测试氢键稳定性
3. **Analysis阶段**：统计氢键并生成PyMOL可视化脚本

## 🚀 快速开始

### 完整流程运行

```bash
# 运行完整的三阶段流程
python scripts/run_hbond_detector.py

# 或者只运行特定阶段
python scripts/run_hbond_detector.py --stage wrapper    # 只运行Wrapper
python scripts/run_hbond_detector.py --stage shaker     # 只运行Shaker
python scripts/run_hbond_detector.py --stage analysis  # 只运行Analysis
```

### 单独分析已有结构

```bash
# 分析已有的PDB文件
python scripts/analyze_hbonds.py data/complex.pdb -o results --csv
```

## 📋 详细使用说明

### 1. 配置文件设置

编辑 `scripts/config.yml` 文件，确保路径正确：

```yaml
paths:
  working_dir: ".."
  autodock4: "autodock4"
  autogrid4: "autogrid4"
  gmx: "gmx"

inputs:
  receptor_pdbqt: "pdbqt/protein.pdbqt"
  ligand_pdbqt: "wrapper/intermediate1_3d.pdbqt"
  ligand_types: "A C N NA OA"

wrapper:
  seeds: [101, 202, 303, 404, 505]
  output_dir: "autodock_runs"

gromacs:
  workdir: "gmx"
```

### 2. Wrapper阶段

生成覆盖蛋白表面的配体群：

```bash
python scripts/run_autodock_batch.py
```

输出：
- `autodock_runs/wrapped_complex.pdb` - 包含所有配体的复合物

### 3. Shaker阶段

通过模拟退火测试氢键稳定性：

```bash
python scripts/washing_cycle.py gmx/
```

关键参数：
- 温度控制：300K → 323K → 300K
- 主链限制：启用 (-DPOSRES)
- 位移阈值：6.0Å

### 4. Analysis阶段

统计氢键并生成可视化：

```bash
python scripts/analyze_hbonds.py gmx/final_complex.pdb -o hbond_results --csv
```

输出：
- `hbond_results.pml` - PyMOL可视化脚本
- `hbond_results.csv` - 详细分析结果

## 📊 结果解读

### WnS评分系统

```
WnS Score = 相互作用能 (E_inter) - (氢键数量 × 2.0)
```

- **分数越低越好**
- 相互作用能：负值表示有利结合
- 氢键贡献：每个氢键减少2.0分

### 氢键判定标准

- **距离**：供体-受体 < 3.5Å
- **角度**：供体-氢-受体 > 120°

### PyMOL可视化

双击 `.pml` 文件或运行：
```bash
pymol hbond_results.pml
```

可视化元素：
- 白色：蛋白质表面
- 青色：配体
- 黄色虚线：氢键

## 🔧 高级配置

### 自定义氢键参数

```bash
python scripts/analyze_hbonds.py complex.pdb -d 3.2 -a 130.0
```

- `-d`：距离阈值（默认3.5Å）
- `-a`：角度阈值（默认120°）

### 批量处理

```bash
# 处理多个PDB文件
for pdb in *.pdb; do
    python scripts/analyze_hbonds.py $pdb -o ${pdb%.pdb}_results --csv
done
```

## 🐛 故障排除

### 常见问题

1. **AutoGrid报错**
   - 检查 `gpf_template.txt` 中的参数文件路径
   - 确保AD4_parameters.dat在正确位置

2. **GROMACS错误**
   - 检查 `annealing.mdp` 中的时间参数
   - 确保拓扑文件正确

3. **氢键分析无结果**
   - 检查配体残基名称（默认为"LIG"）
   - 确认PDB文件格式正确

### 调试模式

```bash
# 干运行（只打印命令不执行）
python scripts/run_hbond_detector.py --dry-run

# 详细输出
python scripts/analyze_hbonds.py complex.pdb -v
```

## 📈 性能优化

### 并行处理

- AutoDock：使用多个seeds并行
- GROMACS：调整 `ntmpi` 和 `ntomp` 参数
- 分析：使用NumPy加速计算

### 内存优化

- 大型系统：分块处理
- 轨迹分析：使用时间窗口

## 📚 参考文献

1. Wrap 'n' Shake算法原理
2. 氢键几何判定标准
3. 模拟退火参数优化

## 🤝 贡献指南

欢迎提交问题报告和功能请求！

## 📄 许可证

本项目遵循MIT许可证。