import os
import shutil
import sys
sys.path.append(r"/proj/c.zihao/work1/function/")
import SWEET

dataset_name = "LUAD"
researchAim = "4noise/15"
base_dir = f"/proj/c.zihao/work1/1NT/{researchAim}"
save_path = f"{base_dir}/SWEET/{dataset_name}"
exprSetFile = f"{base_dir}/data/LUAD_exprSet_noise.csv"
link_file = "/proj/c.zihao/work1/1NT/1data/string/links.csv"

print(f"#========== {dataset_name} ==========")

if os.path.exists(save_path):
    shutil.rmtree(save_path)
os.makedirs(save_path)

SWEET.SWEETcal(
    exprSetFile=exprSetFile,
    link_file=link_file,
    save_path=save_path
)
