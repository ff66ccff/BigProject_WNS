"""
集成示例：展示如何在主流程中使用 StateManager 和 GmxRunner
"""

from state_manager import StateManager
from gmx_runner import run_gmx_mdrun_safe

def run_wrapper():
    """
    Wrapper 循环续传示例
    """
    state = StateManager()
    
    # 1. 读取断点
    start_cycle = state.get("wrapper_cycle", 0)
    max_cycles = 15
    
    print(f"Checking progress... Resuming from Cycle {start_cycle}")

    # 2. 从断点开始循环
    for i in range(start_cycle, max_cycles):
        # ... (运行 AutoGrid, AutoDock, Masking) ...
        
        # 3. 只有本轮全部成功完成后，才保存状态
        # 这样如果中间挂了，下次还是从本轮开始重跑（覆盖本轮坏数据）
        state.update("wrapper_cycle", i + 1)
        state.update("current_receptor", f"cycle_{i+1}/receptor.pdbqt")

    print("Wrapper finished.")

def run_shaker():
    """
    Shaker 阶段续传示例
    """
    state = StateManager()
    
    # 定义有序的阶段列表
    workflow_stages = ["em", "nvt", "annealing_1", "washing_1", "annealing_2", "final_md"]
    
    # 获取上次完成的阶段
    last_finished = state.get("shaker_stage", None)
    
    # 确定开始索引
    start_index = 0
    if last_finished in workflow_stages:
        start_index = workflow_stages.index(last_finished) + 1
        
    # 遍历执行
    for stage in workflow_stages[start_index:]:
        print(f"=== Running Stage: {stage} ===")
        
        if stage == "final_md":
            # 使用方案二的函数，支持 mdrun 内部续传
            run_gmx_mdrun_safe("final_md_production")
        
        elif stage == "washing_1":
            # 运行清洗脚本
            run_washing_script()
            
        # ... 其他阶段逻辑 ...

        # 完成一个阶段，立即存档
        state.update("shaker_stage", stage)

def run_washing_script():
    """
    示例清洗脚本函数
    """
    print("Running washing cycle...")
    # 实际的清洗逻辑