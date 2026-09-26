import os
import shutil
import sys
sys.path.append("/proj/c.zihao/work1/function/benchmark/")
import SSN

dataset_name = "CRC"
base_dir = "/proj/c.zihao/work1/02analysis/03net/2biogrid/out"
save_path = f"{base_dir}/SSN/{dataset_name}"
shutil.rmtree(save_path, ignore_errors=True)
os.makedirs(save_path)
SSN.SSNcal(
    exprSetFile=f"/proj/c.zihao/work1/01data/exprset/{dataset_name}_exprSet_filtered.csv",
    NormalFile=f"/proj/c.zihao/work1/01data/03TCGAnormal/paired_normals/TCGA-{dataset_name}_paired_normal.csv",
    link_file="/proj/c.zihao/work1/01data/04biogrid/biogrid_link.csv",
    organ=None,
    save_path=save_path,
)
