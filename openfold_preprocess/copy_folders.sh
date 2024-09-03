#!/bin/bash

SOURCE_DIR="/data_0/pdb_mmcif"
TARGET_DIR="/data_0/pdb_mmcif/all_cif_content"

# 创建目标目录
mkdir -p "$TARGET_DIR"

# 查找所有文件夹并拷贝到目标目录，同时显示进度
for subdir in "$SOURCE_DIR"/subfolder_*; do
    # 使用rsync拷贝每个子目录的内容到目标目录
    rsync -ah --progress "$subdir"/ "$TARGET_DIR" | pv -l >/dev/null
done
echo "所有内容已拷贝到 $TARGET_DIR"
