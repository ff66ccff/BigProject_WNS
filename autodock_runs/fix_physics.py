import os

def fix_physics_and_config():
    # 1. 修复参数文件：追加 X 原子定义
    param_file = 'AD4_parameters.dat'
    with open(param_file, 'r') as f:
        content = f.read()
    
    if "atom_par X" not in content:
        print(f"Adding X atom definition to {param_file}...")
        with open(param_file, 'a') as f:
            # 保持与其他参数类似的格式
            f.write("\natom_par X     3.60  0.0001  12.000  -0.00110  0.0   0.0  0   -1  -1  3\n")
    else:
        print(f"X atom already defined in {param_file}.")

    # 2. 重新生成 autogrid.gpf (确保引用了 receptor X 但不生成 map X)
    gpf_content = """parameter_file AD4_parameters.dat
npts 60 60 60
gridfld protein.maps.fld
spacing 0.375
receptor protein.pdbqt
gridcenter 0.0 0.0 0.0
receptor_types A C N NA OA X
ligand_types A C N NA OA
map A.map
map C.map
map N.map
map NA.map
map OA.map
elecmap e.map
dsolvmap d.map
"""
    with open('autogrid.gpf', 'w') as f:
        f.write(gpf_content)
    print("Regenerated autogrid.gpf with correct types.")

    # 3. 清理旧 Map 防止干扰
    os.system("del *.map *.fld *.xyz 2>nul")
    print("Cleaned up old map files.")

if __name__ == "__main__":
    fix_physics_and_config()