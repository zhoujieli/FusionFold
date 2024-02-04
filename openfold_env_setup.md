
# 瑶光开发环境数据挂载
1. 瑶光公开数据集：
   * alphafold_dataset
   * uniprot
2. 已经算好的MSA：
   * 自定义存储 访问协议：NFS 协议地址：10.11.12.156 协议路径：/volume1/readonly/common/public/dataset/msa_alignment_output 孙柠他们自己算的MSA
   * 自定义存储 访问协议：NFS 协议地址：10.11.12.156 协议路径：/volume1/1/self-define/yangfei/jyf/preprocess_data/pdb 官方放出来的MSA
3. 标准存储挂载盘：
   * 标准存储 存储卷：protein_one 容器挂载路径：/data_0/
   * 标准存储 存储卷：protein_one 容器挂载路径：/data_1/
   * 标准存储 存储卷：protein_2 容器挂载路径：/data_2/
   * 标准存储 存储卷：protein_3 容器挂载路径：/data_3/
   * 标准存储 存储卷：protein_4 容器挂载路径：/data_4/

# openfold 继续训练环境搭建相关问题踩坑记录
openfold 仓库地址：https://github.com/aqlaboratory/openfold
## 数据准备
1. 数据下载
2. MSA数据库生成
3. MSA查询
## 安装openfold环境
1. 系统环境安装
2. openfold环境安装
   * 直接使用推荐的方法，创建一个新的conda环境：mamba env create -n openfold_env -f environment.yml
   * 或者在已有的环境上进行更新：conda env update --name openfold --file environment.yml
   * 注意environments.yml中可能出现的问题1：
     * 特定版本(5b838a8)的flash-attention无法正常下载。
       * 解决办法：去flash-attention的[commit repo 界面](https://github.com/Dao-AILab/flash-attention/tree/5b838a8bef78186196244a4156ec35bbb58c337d)上，然后直接下载zip文件。注释掉environment.yml中 "- git+https://github.com/Dao-AILab/flash-attention.git@5b838a8" 的部分。
       		1. 在本地解压缩文件夹，发现csrc/falsh_attn路径下，需要cutlass@319389，所以去cutlass的[commit repo 界面](https://github.com/NVIDIA/cutlass/tree/319a389f42b776fae5701afcb943fc03be5b5c25),下载zip文件，放到对应的文件夹中。
       		2. 然后到flash_attention目录下，然后解压缩到本地，从source进行安装： flash_pip install .
	 * 可能报错的问题2: FileNotFoundError: [Errno 2] No such file or directory: '/home/shenyichong/openfold/scripts/cutlass/CHANGELOG.md'
       * 解决方法：将上述cutlass拷贝到openfold对应的文件夹中：/home/shenyichong/openfold/scripts/cutlass ，然后再重新更新conda环境。
4. 电子
