import os

def dlg_to_pdb(dlg_path, output_pdb):
    """从 DLG 文件中提取最终的 DOCKED 构象"""
    with open(dlg_path, 'r') as f:
        lines = f.readlines()
    
    docked_lines = []
    in_docked = False
    found = False
    
    for line in lines:
        if line.startswith("DOCKED: USER"):
            in_docked = True
            found = True
            continue
        if in_docked:
            if line.startswith("DOCKED: ENDMDL"):
                break
            # DLG 中的原子行以 "DOCKED: " 开头，需要去掉这个前缀
            if line.startswith("DOCKED: ATOM") or line.startswith("DOCKED: HETATM"):
                docked_lines.append(line.replace("DOCKED: ", ""))
    
    if found:
        with open(output_pdb, 'w') as f:
            f.write("MODEL 1\n")
            f.writelines(docked_lines)
            f.write("ENDMDL\n")
        print(f"✅ Extracted PDB to: {output_pdb}")
    else:
        print(f"❌ No DOCKED conformation found in {dlg_path}")

# 执行提取
dlg_to_pdb("wrapper_101.dlg", "wrapper_101_docked.pdb")