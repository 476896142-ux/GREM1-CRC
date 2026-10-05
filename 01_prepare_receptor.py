#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Step 1: 受体制备 —— GREM1 (PDB: 5AEJ)
功能：下载晶体结构、去水、加氢、计算电荷、生成 PDBQT
注意：如果 AutoDockTools 已手动准备好受体，此脚本可跳过
"""

import os
import urllib.request

# 目录设置
BASE_DIR = "D:/Bioinfor/colorectal/MR/virtual_screening"
RECEPTOR_DIR = os.path.join(BASE_DIR, "01_receptor")
os.makedirs(RECEPTOR_DIR, exist_ok=True)

# ===== 1. 下载 PDB 结构 =====
pdb_id = "5AEJ"
pdb_url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
pdb_file = os.path.join(RECEPTOR_DIR, f"{pdb_id}.pdb")

print(f"下载 PDB 结构: {pdb_id}")
urllib.request.urlretrieve(pdb_url, pdb_file)
print(f"已保存: {pdb_file}")

# ===== 2. 读取 PDB 并处理（基础版，复杂操作建议用 AutoDockTools GUI）=====
# 注意：完整的受体制备（去水、加氢、电荷计算）建议用 AutoDockTools 手动完成
# 此脚本仅提供基础处理，最终 PDBQT 请以 AutoDockTools 输出为准

print("\n" + "="*60)
print("重要提示：")
print("完整的受体制备建议使用 AutoDockTools (ADFRsuite) 手动完成：")
print("1. 打开 5AEJ.pdb")
print("2. Edit -> Delete Water（去水）")
print("3. Edit -> Hydrogens -> Add（加极性氢）")
print("4. Edit -> Charges -> Compute Gasteiger（计算电荷）")
print("5. File -> Save -> PDBQT 格式")
print("6. 用 Grid Box 工具定义对接盒子（覆盖 DAN 结构域）")
print("="*60)

# ===== 3. 生成对接配置文件模板 =====
config_content = """receptor = 5AEJ_prepared.pdbqt
ligand = ligand.pdbqt
center_x = 0.0
center_y = 0.0
center_z = 0.0
size_x = 20
size_y = 20
size_z = 20
exhaustiveness = 8
num_modes = 9
energy_range = 3
out = output.pdbqt
log = log.txt
"""

config_file = os.path.join(RECEPTOR_DIR, "config_template.txt")
with open(config_file, "w") as f:
    f.write(config_content)
print(f"\n配置模板已生成: {config_file}")
print("请用 AutoDockTools 确定 center_x/y/z 后修改此文件")
