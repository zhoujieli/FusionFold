# 挂载数据
1. 瑶光公开数据集：
   * alphafold_dataset
   * uniprot
2. 已经算好的MSA：
   * 自定义存储 访问协议：NFS 协议地址：10.11.12.156 协议路径：/volume1/readonly/common/public/dataset/msa_alignment_output 孙柠他们自己算的MSA
   * 自定义存储 访问协议：NFS 协议地址：10.11.12.156 协议路径：/volume1/self-define/yangfei/jyf/preprocess_data/pdb 官方放出来的MSA
3. 标准存储挂载盘：
   * 标准存储 存储卷：protein_one 容器挂载路径：/data_0/
   * 标准存储 存储卷：protein_one 容器挂载路径：/data_1/
   * 标准存储 存储卷：protein_2 容器挂载路径：/data_2/
   * 标准存储 存储卷：protein_3 容器挂载路径：/data_3/
   * 标准存储 存储卷：protein_4 容器挂载路径：/data_4/
