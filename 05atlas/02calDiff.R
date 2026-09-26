rm(list = ls())

base_dir <- "/proj/c.zihao/work1/05atlas"
matrix_file <- "/proj/c.zihao/work1/01NT/01string/02benchmark/02merge/merged_matrices/MOSSN/merged_matrix.csv"
meta_file <- "/proj/c.zihao/work1/01data/03TCGA/metadata.csv"

cancers <- c("BLCA", "BRCA", "CRC", "ESCA", "HNSC",
             "KIRC", "LIHC", "LUAD", "LUSC", "PRAD", "STAD")

data <- read.csv(matrix_file, row.names = 1, check.names = FALSE)

meta <- read.csv(meta_file, row.names = 1)
meta <- data.frame(
  sample = meta$Sample,
  cancer = meta$Type,
  patient = substr(meta$Sample, 1, 12),
  type = substr(meta$Sample, 14, 16)
)
meta <- unique(meta[meta$type %in% c("01A", "11A"), ])

atlas_all <- data.frame()
for (cancer_name in cancers) {
  
  meta_one <- meta[meta$cancer == cancer_name, ]
  samples <- intersect(colnames(data), meta_one$sample)
  patients <- unique(substr(samples, 1, 12))

  tumor_col <- paste0(patients, "_01A")
  normal_col <- paste0(patients, "_11A")
  paired <- tumor_col %in% samples & normal_col %in% samples
  tumor_col <- tumor_col[paired]
  normal_col <- normal_col[paired]

  diff <- as.matrix(data[, tumor_col]) - as.matrix(data[, normal_col])

  p_value <- numeric(nrow(diff))
  for (k in seq_len(nrow(diff))) {
    x <- diff[k, ]
    if (length(unique(round(x, 6))) <= 1) {
      p_value[k] <- 1
    } else {
      p_value[k] <- wilcox.test(x, mu = 0, exact = FALSE)$p.value
    }
  }

  atlas <- data.frame(
    link = rownames(diff),
    n_pairs = ncol(diff),
    delta_median = apply(diff, 1, median),
    delta_mean = apply(diff, 1, mean),
    p_value = p_value
  )
  atlas$fdr <- p.adjust(atlas$p_value, method = "BH")

  tau <- quantile(abs(atlas$delta_median), 0.95)
  atlas$direction <- "no_change"
  atlas$direction[atlas$fdr < 0.05 & atlas$delta_median > tau] <- "gained"
  atlas$direction[atlas$fdr < 0.05 & atlas$delta_median < -tau] <- "lost"
  atlas$tau_used <- as.numeric(tau)
  atlas$cancer <- cancer_name

  atlas <- atlas[order(atlas$link, method = "radix"), ]
  atlas_all <- rbind(atlas_all, atlas)

  message(cancer_name, ": total=", nrow(atlas),
          ", gained=", sum(atlas$direction == "gained"),
          ", lost=", sum(atlas$direction == "lost"))
}

rownames(atlas_all) <- NULL
write.csv(atlas_all, file.path(base_dir, "atlas_all.csv"))
