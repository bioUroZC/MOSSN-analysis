rm(list = ls())

base_dir <- "/proj/c.zihao/work1/1NT/9analysis"
setwd(base_dir)

scripts <- c(
  "1calDiff.R",
  "2plotAtlasSummary.R",
  "3recurrentAtlas.R",
  "4module.R",
  "5plotModuleGeneCounts.R",
  "6module_enrichment.R",
  "7plotEdgePrioritization.R",
  "8edgeCaseStudy.R",
  "9edgeSurvival.R"
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
