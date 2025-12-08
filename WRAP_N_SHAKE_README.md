# Wrap 'n' Shake 算法实现

本文档描述了在本项目中实现的 Wrap 'n' Shake 算法，该算法用于多配体对接和分子动力学模拟。

## 算法概述

Wrap 'n' Shake 是一种多配体对接算法，包含两个主要阶段：

1. **Wrapper 阶段**：串行对接，每个配体对接后屏蔽接触的受体原子
2. **Shaker 阶段**：分子动力学模拟，通过模拟退火和清洗循环去除不稳定的配体

## 实现文件

### 核心模块

1. **`mask_pdbqt.py`** - 原子屏蔽函数
   - 计算受体原子与配体原子的距离
   - 将距离 < 3.5 Å 的受体原子类型改为 X，电荷设为 0.000
   - 实现了 Wrap 'n' Shake 的核心屏蔽机制

2. **`wrap_n_shake_docking.py`** - 串行对接循环
   - 实现对接-屏蔽-对接的串行循环
   - 每次对接后自动屏蔽接触的受体原子
   - 替代了原有的批量并行对接

3. **`washing_cycle.py`** - 清洗循环脚本
   - 运行短时MD模拟（模拟退火）
   - 分析配体位移，去除位移 > 6.0 Å 的配体
   - 自动更新拓扑文件和坐标文件

4. **`run_full_wrap_n_shake.py`** - 完整管道脚本
   - 整合 Wrapper 和 Shaker 阶段
   - 提供完整的端到端流程

### 配置文件

1. **`annealing.mdp`** - 模拟退火参数
   - 温度循环：300K → 353K → 300K
   - 位置限制：保持主链限制 (-DPOSRES)
   - 适用于大分子的 80°C 退火温度

2. **`gpf_template.txt`** - AutoGrid 参数模板
   - 添加了 X 原子参数定义
   - Sc2 参数：ε=10^-4 kcal/mol, Rmin≈3.6 Å

## 使用方法

### 1. 配置文件设置

复制 `config.example.yml` 为 `config.yml` 并根据需要调整参数：

```yaml
# Wrap 'n' Shake 特定参数
shaker:
  n_cycles: 5                    # 清洗循环次数
  displacement_cutoff: 6.0       # 位移阈值 (Å)
  cycle_time: 1.0               # 每次循环的MD时间 (ns)
  ligand_resname: "LIG"         # 配体残基名称
```

### 2. 运行完整管道

```bash
# 运行完整的 Wrap 'n' Shake 管道
python run_full_wrap_n_shake.py --config config.yml

# 仅运行 Wrapper 阶段
python run_full_wrap_n_shake.py --config config.yml --skip-shaker

# 仅运行 Shaker 阶段
python run_full_wrap_n_shake.py --config config.yml --skip-wrapper
```

### 3. 单独运行各个阶段

```bash
# 仅运行 Wrapper 阶段
python wrap_n_shake_docking.py --config config.yml

# 运行单个清洗循环
python washing_cycle.py gmx/ -c 6.0 -t 1.0
```

## 算法特点

### Wrapper 阶段特点

1. **串行对接**：不是并行批量对接，而是串行的"对接-屏蔽-对接"循环
2. **原子屏蔽**：每次对接后，将与配体接触的受体原子屏蔽（类型改为X，电荷归零）
3. **排斥力**：X 原子产生短程排斥力，防止后续配体占据相同位置

### Shaker 阶段特点

1. **模拟退火**：使用温度循环（300K → 353K → 300K）测试配体稳定性
2. **位置限制**：保持蛋白质主链限制，仅允许配体移动
3. **清洗循环**：自动去除位移过大的不稳定配体
4. **拓扑更新**：自动更新 GROMACS 拓扑文件，删除被清洗的配体

## 与传统方法的区别

| 特性           | 传统批量对接 | Wrap 'n' Shake   |
| -------------- | ------------ | ---------------- |
| 对接方式       | 并行批量     | 串行循环         |
| 受体处理       | 静态不变     | 动态屏蔽         |
| 配体间相互作用 | 忽略         | 通过屏蔽考虑     |
| MD模拟         | 标准MD       | 模拟退火+清洗    |
| 结果筛选       | 基于对接能   | 基于动力学稳定性 |

## 输出文件

### Wrapper 阶段输出

- `autodock_runs/docked_ligand_*.pdbqt` - 各个对接循环的最佳构象
- `autodock_runs/receptor_masked_*.pdbqt` - 屏蔽后的受体文件

### Shaker 阶段输出

- `gmx/npt.gro` - 最终的复合物结构
- `gmx/topol.top` - 更新后的拓扑文件
- `gmx/anneal.xtc` - MD轨迹文件

## 注意事项

1. **X 原子参数**：确保在 AutoGrid 参数中正确定义 X 原子
2. **位置限制**：Shaker 阶段必须保持蛋白质的位置限制
3. **配体残基名**：确保配体残基名在拓扑文件中正确设置
4. **位移阈值**：可根据系统特性调整位移阈值（默认6.0 Å）

## 参考文献

Wrap 'n' Shake 算法的详细描述请参考原始论文。本实现遵循了论文中描述的核心算法原理，并针对 GROMACS 和 AutoDock 环境进行了优化。