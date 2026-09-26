import os
import shutil
import sys
sys.path.append("/proj/c.zihao/work1/function/benchmark/")
import SWEET

dataset_name = "BLCA"
base_dir = "/proj/c.zihao/work1/02analysis/03net/2intact/out"
save_path = f"{base_dir}/SWEET/{dataset_name}"
shutil.rmtree(save_path, ignore_errors=True)
os.makedirs(save_path)
SWEET.SWEETcal(
    exprSetFile=f"/proj/c.zihao/work1/01data/exprset/{dataset_name}_exprSet_filtered.csv",
    link_file="/proj/c.zihao/work1/01data/07intact/intact_link.csv",
    save_path=save_path,
)
