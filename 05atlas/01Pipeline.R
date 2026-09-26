rm(list = ls())

base_dir <- "/proj/c.zihao/work1/05atlas"
setwd(base_dir)

scripts <- c(
  "02calDiff.R",
  "03recurrent.R",
  "04module.R",
  "05enrichment.R",
  "06plot.R"
)

for (script_name in scripts) {
  message("\n========================================")
  message("Running ", script_name)
  message("========================================")
  script_path <- file.path(base_dir, script_name)
  script_env <- new.env(parent = globalenv())
  source(script_path, echo = FALSE, chdir = TRUE, local = script_env)
}

message("\nAtlas pipeline completed.")
