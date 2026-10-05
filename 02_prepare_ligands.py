#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Step 2: 化合物库准备
功能：从 SMILES 生成 3D 构象，转换为 PDBQT 格式
输入：化合物 SMILES 列表（CSV 格式，含 name 和 smiles 列）
输出：PDBQT 文件
依赖：rdkit, meeko, openbabel
"""

import os
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from meeko import MoleculePreparation, PDBQTWriterLegacy

# 目录设置
BASE_DIR = "D:/Bioinfor/colorectal/MR/virtual_screening"
LIGAND_DIR = os.path.join(BASE_DIR, "02_ligands")
PREPARED_DIR = os.path.join(LIGAND_DIR, "prepared")
os.makedirs(PREPARED_DIR, exist_ok=True)

# ===== 化合物库来源 =====
# 方案 A：DrugBank 已批准药物（推荐先做这个，~2500 个，速度快）
#   下载地址：https://go.drugbank.com/releases/latest（需注册）
#   或从 https://www.drugbank.ca/releases 下载 approved drugs SMILES

# 方案 B：TCMSP 中药活性成分（~13000 个）
#   下载地址：https://tcmsp-e.com/
#   筛选条件：OB >= 30%, DL >= 0.18

# 方案 C：ZINC 类药化合物（~10000 个，可直接下载 PDBQT）
#   下载地址：https://zinc.docking.org/substances/subsets/drug-like/
#   选择 "PDBQT" 格式直接下载，无需转换

# ===== 示例：从 SMILES 列表准备化合物 =====
# 输入文件格式：CSV，两列 name, smiles
# 示例：
# name,smiles
# Aspirin,CC(=O)Oc1ccccc1C(=O)O
# Ibuprofen,CC(C)Cc1ccc(cc1)C(C)C(=O)O

smiles_file = os.path.join(LIGAND_DIR, "ligands_smiles.csv")

if not os.path.exists(smiles_file):
    print(f"未找到 SMILES 文件: {smiles_file}")
    print("请先下载化合物库并整理为 CSV 格式（name,smiles 两列）")
    print("\n推荐下载来源：")
    print("1. DrugBank Approved: https://go.drugbank.com/")
    print("2. TCMSP: https://tcmsp-e.com/")
    print("3. ZINC (直接下载 PDBQT): https://zinc.docking.org/")
    exit(1)

# 读取 SMILES
df = pd.read_csv(smiles_file)
print(f"读取化合物数: {len(df)}")

# ===== 类药性筛选（Lipinski 五规则）=====
def calc_lipinski(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    hbd = Descriptors.NumHDonors(mol)
    hba = Descriptors.NumHAcceptors(mol)
    # Lipinski: MW<=500, LogP<=5, HBD<=5, HBA<=10
    passed = (mw <= 500) and (logp <= 5) and (hbd <= 5) and (hba <= 10)
    return pd.Series([mw, logp, hbd, hba, passed])

df[["MW", "LogP", "HBD", "HBA", "Lipinski"]] = df["smiles"].apply(calc_lipinski)
df_filtered = df[df["Lipinski"] == True].reset_index(drop=True)
print(f"通过 Lipinski 筛选: {len(df_filtered)} / {len(df)}")

# ===== 生成 3D 构象并转换为 PDBQT =====
preparer = MoleculePreparation()
success_count = 0
fail_count = 0

for idx, row in df_filtered.iterrows():
    name = str(row["name"]).replace(" ", "_")
    smiles = row["smiles"]

    try:
        # SMILES -> 3D 构象
        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, randomSeed=42)
        AllChem.MMFFOptimizeMolecule(mol, maxIters=500)

        # RDKit Mol -> PDBQT
        mol_setup = preparer.prepare(mol)
        pdbqt_string = PDBQTWriterLegacy.write_string(mol_setup)[0]

        # 保存
        out_file = os.path.join(PREPARED_DIR, f"{name}.pdbqt")
        with open(out_file, "w") as f:
            f.write(pdbqt_string)

        success_count += 1
        if success_count % 100 == 0:
            print(f"已处理: {success_count} / {len(df_filtered)}")

    except Exception as e:
        fail_count += 1
        continue

print(f"\n转换完成: 成功 {success_count}, 失败 {fail_count}")
print(f"PDBQT 文件保存在: {PREPARED_DIR}")

# 保存筛选后的化合物列表
df_filtered.to_csv(os.path.join(LIGAND_DIR, "ligands_filtered.csv"), index=False)
