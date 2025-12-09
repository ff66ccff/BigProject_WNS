#!/usr/bin/env python3
"""
模拟AutoGrid功能的Python脚本
由于WSL网络问题，我们创建一个临时的grid map文件
"""

import os
import sys
import numpy as np
from struct import pack

def create_dummy_map_files():
    """创建临时的grid map文件"""
    
    # 读取FLD文件获取grid信息
    with open('protein.maps.fld', 'r') as f:
        lines = f.readlines()
    
    # 解析grid信息
    grid_info = {}
    for line in lines:
        if line.startswith('grid'):
            parts = line.strip().split()
            if len(parts) >= 2:
                grid_info[parts[0]] = parts[1]
    
    # 获取grid尺寸
    npts = [60, 60, 60]  # 默认值
    spacing = 0.375
    center = [0.0, 0.0, 0.0]
    
    for line in lines:
        if 'npts' in line:
            npts = list(map(int, line.strip().split()[1:4]))
        elif 'spacing' in line:
            spacing = float(line.strip().split()[1])
        elif 'gridcenter' in line:
            center = list(map(float, line.strip().split()[1:4]))
    
    print(f"Grid info: npts={npts}, spacing={spacing}, center={center}")
    
    # 创建map文件列表
    map_files = ['A.map', 'C.map', 'N.map', 'NA.map', 'OA.map', 'e.map', 'd.map']
    
    # 为每个map文件创建简单的数据
    for map_file in map_files:
        create_map_file(map_file, npts, spacing, center)
    
    print("✅ Created dummy map files")

def create_map_file(filename, npts, spacing, center):
    """创建单个map文件"""
    
    # AutoGrid map文件格式
    # 这是一个简化的实现，创建一个基本的能量场
    
    # 计算grid点总数
    total_points = npts[0] * npts[1] * npts[2]
    
    # 创建简单的能量值（大部分为0，中心区域为负值）
    data = []
    center_x, center_y, center_z = npts[0]//2, npts[1]//2, npts[2]//2
    radius = min(npts) // 4
    
    for z in range(npts[2]):
        for y in range(npts[1]):
            for x in range(npts[0]):
                # 计算到中心的距离
                dist = np.sqrt((x - center_x)**2 + (y - center_y)**2 + (z - center_z)**2)
                
                if dist < radius:
                    # 在中心区域设置负值（有利结合）
                    value = -1.0 * (1.0 - dist/radius)
                else:
                    # 外部区域为0
                    value = 0.0
                
                data.append(value)
    
    # 写入map文件
    with open(filename, 'wb') as f:
        # 写入AutoGrid格式的头部
        f.write(b'GRID\0')
        f.write(pack('i', npts[0]))
        f.write(pack('i', npts[1]))
        f.write(pack('i', npts[2]))
        f.write(pack('d', spacing))
        f.write(pack('d', center[0]))
        f.write(pack('d', center[1]))
        f.write(pack('d', center[2]))
        
        # 写入数据
        for value in data:
            f.write(pack('f', value))

if __name__ == "__main__":
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    create_dummy_map_files()