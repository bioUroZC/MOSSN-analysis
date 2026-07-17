import os
import shutil
import sys

PROJ_ROOT = os.environ.get("MOSSN_ROOT", "/proj/c.zihao/work1")

sys.path.append(PROJ_ROOT + "/function/")
import PPIXpress

dataset_name = "LUAD"
researchAim = "5noise/15"
base_dir = f"{PROJ_ROOT}/1NT/{researchAim}"
save_path = f"{base_dir}/PPIXpress/{dataset_name}"
exprSetFile = f"{base_dir}/data/LUAD_exprSet_noise.csv"
link_file = PROJ_ROOT + "/1NT/1data/string/links.csv"

print(f"#========== {dataset_name} ==========")

if os.path.exists(save_path):
    shutil.rmtree(save_path)
os.makedirs(save_path)

PPIXpress.PPIXpressCal(
    exprSetFile=exprSetFile,
    link_file=link_file,
    save_path=save_path
)
