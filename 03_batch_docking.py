#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Step 3: 批量分子对接（AutoDock Vina）
功能：批量对接所有化合物，提取最佳结合能
支持：多进程并行、断点续跑
依赖：vina（需加入 PATH）、pandas
"""

import os
import subprocess
import pandas as pd
from multiprocessing import Pool
from functools import partial

# ===== 配置 =====
BASE_DIR = "D:/Bioinfor/colorectal/MR/virtual_screening"
RECEPTOR_FILE = os.path.join(BASE_DIR, "01_receptor", "5AEJ_prepared.pdbqt")
LIGAND_DIR = os.path.join(BASE_DIR, "02_ligands", "prepared")
OUTPUT_DIR = os.path.join(BASE_DIR, "03_docking", "output")
LOG_DIR = os.path.join(BASE_DIR, "03_docking", "logs")
RESULT_FILE = os.path.join(BASE_DIR, "04_results", "docking_scores.csv")

# 对接盒子参数（用 AutoDockTools 确定后修改）
CENTER_X = 0.0
CENTER_Y = 0.0
CENTER_Z = 0.0
SIZE_X = 20
SIZE_Y = 20
SIZE_Z = 20
EXHAUSTIVENESS = 8
NUM_MODES = 9

# 并行进程数
N_CORES = 8

os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(LOG_DIR, exist_ok=True)
os.makedirs(os.path.dirname(RESULT_FILE), exist_ok=True)

# ===== 检查受体文件 =====
if not os.path.exists(RECEPTOR_FILE):
    print(f"错误：未找到受体文件 {RECEPTOR_FILE}")
    print("请先用 AutoDockTools 准备受体并保存为 PDBQT")
    exit(1)

# ===== 获取所有配体文件 =====
ligand_files = [f for f in os.listdir(LIGAND_DIR) if f.endswith(".pdbqt")]
print(f"待对接化合物数: {len(ligand_files)}")

# ===== 单化合物对接函数 =====
def dock_single(ligand_file, receptor_file, output_dir, log_dir,
                center_x, center_y, center_z, size_x, size_y, size_z,
                exhaustiveness, num_modes):
    ligand_name = os.path.splitext(ligand_file)[0]
    ligand_path = os.path.join(LIGAND_DIR, ligand_file)
    out_file = os.path.join(output_dir, f"{ligand_name}_out.pdbqt")
    log_file = os.path.join(log_dir, f"{ligand_name}_log.txt")

    # 已对接则跳过（断点续跑）
    if os.path.exists(out_file):
        return None

    # 构建 Vina 命令
    cmd = [
        "vina",
        "--receptor", receptor_file,
        "--ligand", ligand_path,
        "--center_x", str(center_x),
        "--center_y", str(center_y),
        "--center_z", str(center_z),
        "--size_x", str(size_x),
        "--size_y", str(size_y),
        "--size_z", str(size_z),
        "--exhaustiveness", str(exhaustiveness),
        "--num_modes", str(num_modes),
        "--out", out_file,
        "--log", log_file
    ]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)

        # 从日志提取最佳结合能
        best_score = None
        if os.path.exists(log_file):
            with open(log_file, "r") as f:
                lines = f.readlines()
                for line in lines:
                    if line.startswith("   1 "):
                        parts = line.split()
                        if len(parts) >= 2:
                            best_score = float(parts[1])
                        break

        return {
            "ligand": ligand_name,
            "best_score": best_score,
            "status": "success" if best_score else "no_score"
        }
    except subprocess.TimeoutExpired:
        return {"ligand": ligand_name, "best_score": None, "status": "timeout"}
    except Exception as e:
        return {"ligand": ligand_name, "best_score": None, "status": f"error: {str(e)}"}

# ===== 批量对接（并行）=====
print(f"\n开始批量对接，使用 {N_CORES} 个进程...")
print(f"对接盒子: center=({CENTER_X},{CENTER_Y},{CENTER_Z}), size=({SIZE_X},{SIZE_Y},{SIZE_Z})")

dock_func = partial(
    dock_single,
    receptor_file=RECEPTOR_FILE,
    output_dir=OUTPUT_DIR,
    log_dir=LOG_DIR,
    center_x=CENTER_X, center_y=CENTER_Y, center_z=CENTER_Z,
    size_x=SIZE_X, size_y=SIZE_Y, size_z=SIZE_Z,
    exhaustiveness=EXHAUSTIVENESS, num_modes=NUM_MODES
)

with Pool(processes=N_CORES) as pool:
    results = pool.map(dock_func, ligand_files)

# ===== 汇总结果 =====
results = [r for r in results if r is not None]
df_results = pd.DataFrame(results)
df_results = df_results.sort_values("best_score", na_position="last")
df_results.to_csv(RESULT_FILE, index=False)

print(f"\n对接完成！")
print(f"成功: {len(df_results[df_results['status']=='success'])}")
print(f"失败: {len(df_results[df_results['status']!='success'])}")
print(f"结果保存在: {RESULT_FILE}")
print(f"\nTOP 10 化合物:")
print(df_results.head(10).to_string(index=False))
