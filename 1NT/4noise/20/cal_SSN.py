import os
import shutil
import sys
sys.path.append(r"/proj/c.zihao/work1/function/")
import SSN

dataset_name = "LUAD"
researchAim = "4noise/20"
base_dir = f"/proj/c.zihao/work1/1NT/{researchAim}"
save_path = f"{base_dir}/SSN/{dataset_name}"

NormalFile = f"{base_dir}/data/Lung_normal_reference.csv"
exprSetFile = f"{base_dir}/data/LUAD_exprSet_noise.csv"
link_file = "/proj/c.zihao/work1/1NT/1data/string/links.csv"

print(f"#========== {dataset_name} ==========")

if os.path.exists(save_path):
    shutil.rmtree(save_path)
os.makedirs(save_path)

SSN.SSNcal(
    exprSetFile=exprSetFile,
    NormalFile=NormalFile,
    link_file=link_file,
    organ="Lung",
    save_path=save_path
)
