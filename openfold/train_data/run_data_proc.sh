#!/bin/bash

# 定义目录变量
pdb_dir="./pdb_dir"            # 输入：指定pdb文件的目录
mmcif_dir="./mmcif_dir"        # 输出：指定输出mmcif文件的目录
mmcif_cache_path="./mmcif_cache.json"    # 输出：mmcif缓存文件的路径
chain_data_cache_path="./chain_data_cache.json" # 输出：chain数据缓存文件的路径

# 第一步：转换pdb到mmcif
echo "Step 1: Converting PDB to mmCIF..."
python /data_1/data_process/pdb2cif_3.py $pdb_dir $mmcif_dir --no_workers 16

# 第二步：生成mmcif缓存文件
echo "Step 2: Generating mmCIF cache..."
python /home/fanminzhi/generate_mmcif_cache.py $pdb_dir $mmcif_cache_path --no_workers 16

# 第三步：生成chain数据缓存文件
echo "Step 3: Generating chain data cache..."
python /data_1/openfold/scripts/generate_chain_data_cache.py $pdb_dir/ $chain_data_cache_path --no_workers 16

echo "Data processing completed!"
