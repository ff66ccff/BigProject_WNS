# 断点续传系统使用指南

本项目已集成完整的断点续传系统，支持在任何阶段中断后恢复执行。

## 核心组件

### 1. StateManager (utils/state_manager.py)
- 持久化保存工作流进度到 JSON 文件
- 支持读取和更新状态
- 自动处理文件损坏情况

### 2. GmxRunner (utils/gmx_runner.py)
- 自动检测检查点文件(.cpt)
- 支持从断点续传或新建运行
- 错误处理和状态提示

## 使用方法

### 1. 完整流程 (run_full_pipeline_checkpoint.py)
```bash
# 运行完整流程
python scripts/run_full_pipeline_checkpoint.py --config config.yml

# 从中断点继续运行（自动检测）
python scripts/run_full_pipeline_checkpoint.py --config config.yml

# 重置所有进度，从头开始
python scripts/run_full_pipeline_checkpoint.py --config config.yml --reset
```

### 2. Wrap 'n' Shake 流程 (run_full_wrap_n_shake.py)
```bash
# 运行 Wrap 'n' Shake 流程
python scripts/run_full_wrap_n_shake.py --config config.yml

# 从中断点继续运行
python scripts/run_full_wrap_n_shake.py --config config.yml

# 重置进度
python scripts/run_full_wrap_n_shake.py --config config.yml --reset
```

### 3. GROMACS 流程 (gromacs_pipeline.py)
```bash
# 运行 GROMACS 流程
python scripts/gromacs_pipeline.py complex.pdb gmx_output ligand.mol2

# 自动从检查点恢复（无需额外参数）
python scripts/gromacs_pipeline.py complex.pdb gmx_output ligand.mol2
```

## 状态文件

系统会在工作目录中创建以下状态文件：
- `workflow_state.json` - 主工作流状态
- `gromacs_state.json` - GROMACS 流程状态
- `pipeline_state.json` - 完整流程状态

## 断点续传特性

### 流程级续传
- 跳过已完成的 Wrapper 循环
- 跳过已完成的 Shaker 阶段
- 自动检测并从下一阶段开始

### 计算级续传
- GROMACS 模拟中断后基于 .cpt 文件续传
- 自动添加 `-cpi` 和 `-append` 参数
- 避免重新运行已完成的部分

## 示例场景

### 场景 1：电源中断
```bash
# 电源中断后，直接重新运行命令
python scripts/run_full_pipeline_checkpoint.py --config config.yml
# 系统会自动检测进度并从中断点继续
```

### 场景 2：手动暂停
```bash
# 按 Ctrl+C 暂停后，可以安全地重新运行
python scripts/run_full_wrap_n_shake.py --config config.yml
# 系统会从当前阶段继续
```

### 场景 3：GROMACS 模拟中断
```bash
# GROMACS 模拟被系统杀死后
python scripts/gromacs_pipeline.py complex.pdb gmx_output ligand.mol2
# 会自动从 .cpt 文件恢复
```

## 注意事项

1. **不要手动编辑状态文件**，除非您清楚了解其结构
2. **使用 --reset 参数**会清除所有进度，请谨慎使用
3. **确保工作目录不变**，移动工作目录会导致状态丢失
4. **配置文件修改**后，建议使用 --reset 重新开始

## 故障排除

### 问题：状态文件损坏
```bash
# 解决方案：使用 --reset 重置
python scripts/run_full_pipeline_checkpoint.py --config config.yml --reset
```

### 问题：从错误的阶段开始
```bash
# 解决方案：手动编辑状态文件或使用 --reset
python scripts/run_full_pipeline_checkpoint.py --config config.yml --reset
```

### 问题：GROMACS 检查点文件丢失
```bash
# 解决方案：删除 gmx_state.json 重新开始
rm gmx_state.json
python scripts/gromacs_pipeline.py complex.pdb gmx_output ligand.mol2
```