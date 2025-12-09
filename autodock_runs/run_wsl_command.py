#!/usr/bin/env python3
"""
通过Python运行WSL命令，绕过代理问题
"""
import subprocess
import os
import sys

def run_wsl_command(command):
    """在WSL中运行命令，禁用代理"""
    wsl_command = [
        'wsl', '-d', 'Ubuntu', 'bash', '-c', 
        f'unset http_proxy https_proxy all_proxy && {command}'
    ]
    
    try:
        result = subprocess.run(
            wsl_command,
            capture_output=True,
            text=True,
            timeout=300  # 5分钟超时
        )
        
        print(f"Exit code: {result.returncode}")
        if result.stdout:
            print(f"STDOUT:\n{result.stdout}")
        if result.stderr:
            print(f"STDERR:\n{result.stderr}")
            
        return result.returncode == 0
        
    except subprocess.TimeoutExpired:
        print("Command timed out")
        return False
    except Exception as e:
        print(f"Error running command: {e}")
        return False

if __name__ == "__main__":
    # 切换到当前目录
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    
    # 运行autogrid4
    print("Running AutoGrid...")
    success = run_wsl_command(
        "cd /mnt/d/PersonalFile/Documents/BigProject/SRC/project/autodock_runs && /usr/bin/autogrid4 -p autogrid.gpf -l autogrid.log"
    )
    
    if success:
        print("✅ AutoGrid completed successfully")
        
        # 运行autodock4
        print("\nRunning AutoDock...")
        success = run_wsl_command(
            "cd /mnt/d/PersonalFile/Documents/BigProject/SRC/project/autodock_runs && /usr/bin/autodock4 -p wrapper_101.dpf -l wrapper_101.dlg"
        )
        
        if success:
            print("✅ AutoDock completed successfully")
        else:
            print("❌ AutoDock failed")
    else:
        print("❌ AutoGrid failed")