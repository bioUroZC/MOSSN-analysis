rm(list = ls())

PROJ_ROOT <- Sys.getenv("MOSSN_ROOT", "/proj/c.zihao/work1")


base_dir <- paste0(PROJ_ROOT, "/1NT/9analysis")
setwd(base_dir)

scripts <- c(
  "1calDiff.R",
  "2recurrentAtlas.R",
  "3module.R",
  "4moduleEnrichment.R",
  "5plotPerturbationAtlas.R"
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
