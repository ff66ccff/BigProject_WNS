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

### Step 5: AutoGrid/AutoDock 批量对接
- 脚本：`run_autodock_batch.py`
- 功能：生成网格、执行多 seed 对接
- 输出：`autodock_runs/wrapper_<seed>.dlg`

### Step 6: 构建复合体
- 脚本：`build_complex.py`
- 功能：合并蛋白与最佳配体 pose，过滤原子冲突
- 输出：`complex/complex_filtered.pdb`

### Step 7: GROMACS 分子动力学
- 脚本：`gromacs_pipeline.sh`（在 WSL 中执行）
- 流程：pdb2gmx → solvate → genion → em → nvt → npt → md
- 输出：`gmx/` 目录下的轨迹和拓扑文件

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

| 脚本                    | 用途              | 运行环境      |
| ----------------------- | ----------------- | ------------- |
| `run_full_pipeline.py`  | 一键自动化        | Windows       |
| `preprocess_pdb.py`     | PDB 清洗          | Windows       |
| `run_autodock_batch.py` | AutoDock 批量对接 | Windows → WSL |
| `build_complex.py`      | 复合体构建        | Windows       |
| `ligand_param.sh`       | 配体参数化        | WSL           |
| `gromacs_pipeline.sh`   | MD 模拟           | WSL           |
| `post_md_analysis.py`   | MD 后处理         | Windows/WSL   |
| `cluster.sh`            | 轨迹聚类          | WSL           |
| `update_manifest.py`    | 更新 manifest     | Windows       |
| `package_results.py`    | 归档打包          | Windows       |
