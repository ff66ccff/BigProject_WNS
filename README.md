# 项目使用指南

本目录实现了基于论文 "Sc2" 方法的 AutoDock 对接与 GROMACS 分子动力学模拟的完整自动化流程。建议阅读下方说明后再运行脚本，以确保与本地环境一致。

---

## 目录结构说明

```
project/
├── data/           # 原始输入文件（蛋白 PDB、配体 PDB/mol2 等）
├── pdb/            # 经 preprocess_pdb.py 清洗后的蛋白结构
├── pdbqt/          # 由 MGLTools 生成的受体 PDBQT 文件
├── wrapper/        # 配体 PDBQT 及 AutoDock 包裹相关文件
├── autodock_runs/  # AutoGrid/AutoDock 的日志与 dlg/pdbqt 输出
├── complex/        # 蛋白-配体复合体及几何筛选后的结果
├── gmx/            # GROMACS 构建体系的中间文件
├── md/             # 生产阶段 MD 轨迹、能量文件等
├── mdp/            # GROMACS 参数文件（em/nvt/npt/md.mdp）
├── analysis/       # MD 后处理分析结果
├── manifest/       # 运行记录（run-manifest.yml）与可复现性打包
└── scripts/        # 全部自动化脚本、模板与配置
    ├── config.yml           # 用户配置（从 config.example.yml 复制）
    ├── config.example.yml   # 配置模板
    ├── run_full_pipeline.py # 一键自动化脚本
    ├── preprocess_pdb.py    # PDB 清洗
    ├── prepare_pdbqt.ps1    # MGLTools 受体转换（PowerShell）
    ├── run_autodock_batch.py# AutoGrid/AutoDock 批量对接
    ├── ligand_param.sh      # AmberTools 配体参数化（WSL）
    ├── build_complex.py     # 合并蛋白与配体 pose
    ├── gromacs_pipeline.sh  # GROMACS 全流程（WSL）
    ├── post_md_analysis.py  # MD 后处理（MDAnalysis）
    ├── cluster.sh           # 轨迹聚类
    ├── package_results.py   # 可复现性归档
    ├── update_manifest.py   # 更新 manifest
    └── templates/           # AutoGrid/AutoDock 模板文件
```

---

## 环境与前置要求

| 组件                  | 位置                  | 说明                                     |
| --------------------- | --------------------- | ---------------------------------------- |
| MGLTools              | Windows `D:\MGLTools` | 含 `python.exe` 及 AutoDockTools 脚本    |
| AutoDock4 / AutoGrid4 | WSL2 Ubuntu           | 需在 PATH 中可直接调用                   |
| AmberTools / ACPYPE   | WSL2 Ubuntu           | 配体参数化（可选，若已有力场参数可跳过） |
| GROMACS               | WSL2 Ubuntu           | 需在 PATH 中可直接调用 `gmx`             |
| Python 3              | Windows               | 用于执行项目内 Python 脚本               |
| PyYAML                | Windows Python        | `pip install pyyaml`                     |
| MDAnalysis            | Windows/WSL Python    | `pip install MDAnalysis`（后处理分析用） |

---

## 快速开始

### 1. 准备配置文件

```powershell
cd D:\PersonalFile\Documents\BigProject\SRC\project\scripts
copy config.example.yml config.yml
# 编辑 config.yml，填写实际路径和文件名
```

### 2. 放置输入文件

将原始蛋白 PDB 和配体文件放入 `project/data/` 目录，并在 `config.yml` 中指定文件名：

```yaml
inputs:
  raw_protein_pdb: "data/your_protein.pdb"
  raw_ligand_mol2: "data/your_ligand.pdb"  # 支持 pdb/mol2 格式
```

### 3. 检查配置（Dry-Run）

```powershell
D:\Python\python.exe scripts\run_full_pipeline.py --dry-run
```

确认输出的命令路径正确后再正式运行。

### 4. 正式运行

```powershell
D:\Python\python.exe scripts\run_full_pipeline.py
```

---

## 流程步骤详解

### Step 1: 蛋白 PDB 预处理
- 脚本：`preprocess_pdb.py`
- 功能：标准化链名、去除重复原子记录
- 输出：`pdb/protein_clean.pdb`

### Step 2: 生成受体 PDBQT
- 工具：MGLTools `prepare_receptor4.py`
- 功能：添加 Gasteiger 电荷、生成 PDBQT
- 输出：`pdbqt/protein.pdbqt`

### Step 3: 生成配体 PDBQT
- 工具：MGLTools `prepare_ligand4.py`
- 功能：配体格式转换、添加电荷
- 输出：`wrapper/<ligand>.pdbqt`

### Step 4: 配体参数化（可选）
- 脚本：`ligand_param.sh`（在 WSL 中执行）
- 工具：AmberTools antechamber + ACPYPE
- 功能：生成 GAFF 力场参数
- 注意：若 AmberTools 未安装，此步会报错但不影响对接
 - 配置与可选行为：
   - `config.yml` 中 `ambertools.skip`（布尔，默认 `false`）可以设置为 `true` 来跳过此步骤。
   - 如果你希望在缺少 `antechamber` 时终止整个流程，可以将 `ambertools.fail_on_missing` 设置为 `true`（默认 `false`）。
   - 脚本在执行前会检测 WSL 中是否存在 `antechamber`，若不存在且 `skip` 为 `false`，默认会打印警告并跳过该步骤（不终止流程）。

### Step 5: AutoGrid/AutoDock 批量对接
- 脚本：`run_autodock_batch.py`
- 功能：生成网格、执行多 seed 对接
- 输出：`autodock_runs/wrapper_<seed>.dlg`

### Step 6: 构建复合体
- 脚本：`build_complex.py`
- 功能：合并蛋白与最佳配体 pose，过滤原子冲突
- 输出：`complex/complex_filtered.pdb`

### Step 7: GROMACS 分子动力学
- 脚本：`gromacs_full_auto.sh` v2.1（在 WSL 中执行）**【全自动化脚本 - 已修复】**
- 流程：提取蛋白/配体 → ACPYPE 参数化 → pdb2gmx → 组装 → 盒子 → 溶剂化 → 离子 → 拓扑验证 → EM → NVT → NPT → MD
- 特点：
  - **完全自动化**，无需手动输入任何参数
  - 自动检测已有的配体参数文件（跳过重复计算）
  - 自动处理力场选择、水模型、温度耦合组等
  - 自动创建索引组 (`Protein_UNL`, `Water_and_ions`)
  - **v2.1 修复**：
    - 🔧 修复重复运行导致的拓扑文件累积问题
    - 🔧 在溶剂化前自动重建干净的 `[ molecules ]` 部分
    - 🔧 离子添加后自动验证并修复拓扑一致性
    - 🔧 支持 `Protein_UNL`/`Protein_LIG` 等多种组名
    - 🔧 添加自动清理旧中间文件功能
  - 彩色日志输出，便于追踪进度
- 技术路线对应：
  - Amber ff + TIP3P（根据技术路线第5节）
  - Steepest Descent, Fmax < 1000（根据技术路线第6节）
  - v-rescale + Parrinello-Rahman（根据技术路线第6节）
  - 轨迹采样 2 ps（根据技术路线第6节）
- 输出：`gmx/` 目录下的轨迹和拓扑文件
- 使用方法：
  ```bash
  cd /mnt/d/PersonalFile/Documents/BigProject/SRC/project
  bash scripts/gromacs_full_auto.sh complex/complex_filtered.pdb data/intermediate1.mol2 gmx mdp
  ```
- 环境变量：
  - `CLEAN_OLD=0`：禁用自动清理旧文件（默认启用）
  - `NTMPI=1`：MPI 线程数
  - `NTOMP=8`：OpenMP 线程数

### Step 8: MD 后处理分析
- 脚本：`post_analysis.sh`（在 WSL 中执行）**【新增】**
- 功能（根据技术路线第7节）：
  - 配体占有率密度图 (Occupancy Map)
  - 接触频率统计 (Contact Frequency)
  - 停留时间 (Residence Time)
  - 构象簇分析 (Clustering)
  - RMSD/RMSF 分析
  - 氢键分析
- 使用方法：
  ```bash
  bash scripts/post_analysis.sh gmx
  ```

---

## 配置文件说明 (config.yml)

```yaml
paths:
  working_dir: ".."              # 相对于 scripts/ 的项目根目录
  python: "python"               # Windows Python 可执行文件
  mgltools_root: "D:/MGLTools"   # MGLTools 安装目录
  autodock4: "autodock4"         # WSL 中 autodock4 命令
  autogrid4: "autogrid4"         # WSL 中 autogrid4 命令
  gmx: "gmx"                     # WSL 中 gmx 命令

inputs:
  raw_protein_pdb: "data/MAB_1134c.pdb"
  raw_ligand_mol2: "data/MAB_2301.pdb"
  ligand_types: "A C HD N OA"    # 配体原子类型

autogrid:
  npts: [60, 60, 60]             # 网格点数
  center: [0.0, 0.0, 0.0]        # 网格中心（需根据蛋白调整）
  spacing: 0.375                 # 网格间距 (Å)

wrapper:
  seeds: [101, 202, 303, 404, 505]  # AutoDock 随机种子
```

---

## 常见问题

### Q: MGLTools 找不到 python.exe？
A: 检查 `config.yml` 中 `mgltools_root` 路径，确保该目录下存在 `python.exe` 或 `mglpython.exe`。

### Q: WSL 中 antechamber/autodock4 找不到？
A: 确保已在 WSL 中安装并配置好 PATH。可在 WSL 终端测试：
```bash
which autodock4
which antechamber
```

### Q: 如何跳过配体参数化步骤？
A: 若 AmberTools 未安装，该步骤会报错退出。可在 `config.yml` 中添加 `skip_ambertools: true`（需修改脚本支持），或手动注释掉相关代码。

### Q: AutoGrid 网格中心如何确定？
A: 通常设为蛋白活性位点的几何中心。可用 PyMOL 或 VMD 查看蛋白结构后手动填写。

---

## 可复现性归档

执行以下命令生成归档包：

```powershell
D:\Python\python.exe scripts\package_results.py
```

归档包包含：manifest、脚本、分析结果，方便论文附录或数据共享。

---

## 脚本速查表

| 脚本                     | 用途                 | 运行环境      |
| ------------------------ | -------------------- | ------------- |
| `run_full_pipeline.py`   | 一键自动化           | Windows       |
| `preprocess_pdb.py`      | PDB 清洗             | Windows       |
| `run_autodock_batch.py`  | AutoDock 批量对接    | Windows → WSL |
| `build_complex.py`       | 复合体构建           | Windows       |
| `ligand_param.sh`        | 配体参数化           | WSL           |
| `gromacs_pipeline.sh`    | MD 模拟（旧版）      | WSL           |
| `gromacs_full_auto.sh`   | **MD 全自动化 v2.0** | WSL           |
| `continue_simulation.sh` | 从中断点继续 MD      | WSL           |
| `post_analysis.sh`       | **MD 后处理分析**    | WSL           |
| `post_md_analysis.py`    | MD 后处理 (Python)   | Windows/WSL   |
| `cluster.sh`             | 轨迹聚类             | WSL           |
| `update_manifest.py`     | 更新 manifest        | Windows       |
| `package_results.py`     | 归档打包             | Windows       |

---

## 技术路线对应

本项目严格按照技术路线文档实现：

| 技术路线步骤      | 脚本/工具                   | 关键参数                |
| ----------------- | --------------------------- | ----------------------- |
| 1. 蛋白质结构准备 | `preprocess_pdb.py`         | 格式标准化              |
| 2. PDB→PDBQT 转换 | MGLTools                    | Gasteiger 电荷          |
| 3. Wrapper (Sc2)  | AutoDock                    | ε=1×10⁻⁴, R≈3.6Å        |
| 4. 复合体构建     | `build_complex.py` + ACPYPE | GAFF2 力场              |
| 5. 体系构建       | GROMACS                     | Amber ff, TIP3P, 1.0nm  |
| 6. MD 模拟        | GROMACS                     | v-rescale, P-R, 2ps采样 |
| 7. 后处理分析     | `post_analysis.sh`          | 占有率/接触/聚类        |

### MDP 参数配置

| 参数       | 设置值            | 说明          |
| ---------- | ----------------- | ------------- |
| 力场       | amber99sb-ildn    | Amber ff 系列 |
| 水模型     | TIP3P             | 与 Amber 兼容 |
| 静电       | PME               | 长程静电处理  |
| 能量最小化 | Steepest Descent  | Fmax < 1000   |
| 温控器     | v-rescale         | 300 K         |
| 压控器     | Parrinello-Rahman | 1 bar         |
| 轨迹采样   | 1000 步 = 2 ps    | 根据技术路线  |
