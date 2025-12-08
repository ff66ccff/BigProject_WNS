# 项目使用指南

本目录实现了基于论文 "Wrap 'n' Shake (WnS)" 方法的完整自动化流程，包括重构后的Wrapper覆盖策略、Shaker模拟协议、Washing过滤机制和幸存者评估指标。建议阅读下方说明后再运行脚本，以确保与本地环境一致。

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

### Step 5: Wrapper（重构后的串行迭代屏蔽对接）
- 脚本：`run_autodock_batch.py`（已重构为串行迭代模式）
- 功能：实现论文要求的"单分子层"覆盖策略
- 核心修改：
  - **串行迭代对接**：从"并行独立对接"改为"串行迭代屏蔽对接"
  - **屏蔽机制**：识别已对接配体周围3.5Å范围内的蛋白原子，将其原子类型强制改为"X"
  - **原子"X"参数**：电荷q=0，Lennard-Jones参数ε=10⁻⁴ kcal/mol, Rmin≈3.6Å
  - **终止条件**：当蛋白表面覆盖率达到饱和（剩余空白表面积<1%）时停止
- 输出：`autodock_runs/wrapper_<seed>.dlg`，`wrapper/monolayer.pdbqt`

### Step 6: 构建复合体
- 脚本：`build_complex.py`
- 功能：合并蛋白与最佳配体 pose，过滤原子冲突
- 输出：`complex/complex_filtered.pdb`

### Step 7: Shaker（重构后的受限模拟退火协议）
- 脚本：`gromacs_full_auto.sh` v3.0（在 WSL 中执行）**【重构为WnS协议】**
- 流程：提取蛋白/配体 → ACPYPE 参数化 → pdb2gmx → 组装 → 盒子 → 溶剂化 → 离子 → 拓扑验证 → EM → **MDB（初筛MD）** → **MDBSA（模拟退火MD）** → **Washing（清洗）** → **MDF（全柔性精修）**
- 核心修改：
  - **初筛MD (MDB)**：短时间（5ns）MD，对蛋白主链施加位置限制，仅允许侧链和配体移动
  - **模拟退火MD (MDBSA)**：
    - 保留蛋白主链的位置限制（防止高温下蛋白解折叠）
    - 锯齿状温度变化协议：
      - 小分子 (MW ≤ 300)：升温至 323 K (50°C)
      - 大分子 (MW ≥ 300)：升温至 353 K (80°C)
    - 利用高温加速弱结合配体的解离
  - **Washing（清洗）**：
    - 在每一轮Shaker模拟结束后分析轨迹
    - 判定标准：配体质心移动超过6Å或相互作用能变弱
    - 清洗动作：从系统中删除已解离的配体
    - 循环直至约75%的配体被移除
  - **全柔性精修 (MDF)**：对幸存配体解除限制，进行无约束、300K的常规MD精修
- 输出：`gmx/` 目录下的轨迹和拓扑文件，`md/survivors.pdb`
- 使用方法：
  ```bash
  cd /mnt/d/PersonalFile/Documents/BigProject/SRC/project
  bash scripts/gromacs_full_auto.sh wrapper/monolayer.pdbqt data/ligand.mol2 gmx mdp
  ```
- 环境变量：
  - `CLEAN_OLD=0`：禁用自动清理旧文件（默认启用）
  - `NTMPI=1`：MPI 线程数
  - `NTOMP=8`：OpenMP 线程数
  - `SHAKER_TEMP=323`：小分子温度（默认323K）
  - `SHAKER_TEMP_LARGE=353`：大分子温度（默认353K）
  - `WASHING_THRESHOLD=6.0`：清洗距离阈值（默认6.0Å）
  - `SURVIVAL_RATE=0.25`：幸存率阈值（默认25%）

### Step 8: 幸存者分析与结合位点评估
- 脚本：`post_analysis.sh`（在 WSL 中执行）**【重构为幸存者评估】**
- 功能（根据WnS方法第7节）：
  - **Shaker Rate (SR)**：计算公式 SR = N_initial/N_final，SR值越高说明非特异性位点去除得越干净
  - **幸存者排序**：根据最终幸存配体与蛋白的相互作用能（Interaction Energy）进行排序
  - **聚类分析**：对幸存下来的配体位置进行RMSD聚类，簇中心即为预测的结合位点
  - **结合位点可视化**：生成幸存配体的占有率密度图和结合位点预测
  - **能量分析**：计算幸存配体的结合自由能和相互作用能分解
- 使用方法：
  ```bash
  bash scripts/post_analysis.sh gmx
  ```
- 输出：`analysis/` 目录下的幸存者分析结果、结合位点预测和可视化图表

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

## 修改后的技术路线图

本项目实现了重构后的Wrap 'n' Shake (WnS)技术路线：

| 技术路线步骤       | 脚本/工具                   | 关键参数/修改点               |
| ------------------ | --------------------------- | ----------------------------- |
| 1. 蛋白质结构准备  | `preprocess_pdb.py`         | 格式标准化                    |
| 2. PDB→PDBQT 转换  | MGLTools                    | Gasteiger 电荷                |
| 3. Wrapper（重构） | `run_autodock_batch.py`     | 串行迭代屏蔽对接，原子"X"参数 |
| 4. 复合体构建      | `build_complex.py` + ACPYPE | GAFF2 力场                    |
| 5. 体系构建        | GROMACS                     | Amber ff, TIP3P, 1.0nm        |
| 6. Shaker（重构）  | `gromacs_full_auto.sh`      | 受限模拟退火+Washing机制      |
| 7. 幸存者分析      | `post_analysis.sh`          | Shaker Rate+幸存者聚类        |

### 核心修改点总结

**修改点一：重构 Wrapper 覆盖策略**
- 从"并行独立对接"改为"串行迭代屏蔽对接"
- 实现蛋白原子类型"X"的动态修改和屏蔽机制
- 确保形成真正的"单分子层"覆盖

**修改点二：重构 Shaker 模拟协议**
- 从"常规全柔性MD"改为"受限模拟退火"
- 实现三阶段MD：初筛MD→模拟退火MD→全柔性精修
- 根据分子量设定不同的温度协议

**修改点三：增加 Washing 过滤机制**
- 在Shaker模拟过程中实时分析配体状态
- 自动删除解离的配体，防止干扰后续模拟
- 实现约75%配体移除的终止条件

**修改点四：调整评估指标**
- 从"占有率密度图"改为"幸存者"评估
- 引入Shaker Rate (SR)指标量化非特异性位点去除效果
- 基于幸存配体的聚类预测结合位点

### MDP 参数配置（重构后的WnS协议）

| 参数           | 设置值              | 说明          |
| -------------- | ------------------- | ------------- |
| 力场           | amber99sb-ildn      | Amber ff 系列 |
| 水模型         | TIP3P               | 与 Amber 兼容 |
| 静电           | PME                 | 长程静电处理  |
| 能量最小化     | Steepest Descent    | Fmax < 1000   |
| 初筛MD温控器   | v-rescale           | 300 K         |
| 模拟退火温控器 | v-rescale           | 323-353 K     |
| 压控器         | Parrinello-Rahman   | 1 bar         |
| 轨迹采样       | 1000 步 = 2 ps      | 根据技术路线  |
| 位置限制       | Position Restraints | 蛋白主链限制  |

### 新增配置参数（config.yml）

```yaml
shaker:
  # 模拟退火参数
  small_mol_temp: 323      # 小分子温度 (K)
  large_mol_temp: 353      # 大分子温度 (K)
  mol_weight_threshold: 300 # 分子量阈值 (Da)
  
  # Washing参数
  washing_threshold: 6.0   # 配体移动距离阈值 (Å)
  survival_rate: 0.25      # 幸存率阈值
  
  # MD时间设置
  pre_md_time: 5           # 初筛MD时间 (ns)
  annealing_cycles: 5      # 模拟退火循环数
  final_md_time: 10        # 最终精修MD时间 (ns)

wrapper:
  # 屏蔽对接参数
  masking_radius: 3.5      # 屏蔽半径 (Å)
  x_atom_epsilon: 1e-4     # X原子LJ参数ε
  x_atom_rmin: 3.6         # X原子LJ参数Rmin (Å)
  coverage_threshold: 0.99 # 覆盖率阈值
```
