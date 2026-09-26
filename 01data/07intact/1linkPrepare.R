rm(list = ls())

library(data.table)

data_dir <- "/proj/c.zihao/work1/01data/07intact"
zip_path <- file.path(data_dir, "intact.zip")
output_path <- file.path(data_dir, "intact_link.csv")

normalize_intact_names <- function(dt) {
  setnames(dt, old = names(dt), new = sub("^#", "", names(dt)))
  dt
}

extract_first_match <- function(x, pattern) {
  hit <- regexpr(pattern, x, perl = TRUE)
  out <- rep(NA_character_, length(x))
  ok <- hit > 0
  out[ok] <- regmatches(x, hit)
  out
}

if (!file.exists(zip_path)) stop("IntAct archive not found: ", zip_path)
members <- unzip(zip_path, list = TRUE)$Name
if (!("intact.txt" %in% members)) stop("intact.txt is absent from: ", zip_path)

intact_dt <- fread(
  cmd = sprintf("unzip -p %s %s", shQuote(zip_path), shQuote("intact.txt")),
  sep = "\t", header = TRUE, quote = "", fill = TRUE, showProgress = TRUE
)
normalize_intact_names(intact_dt)

edge_dt <- intact_dt[, .(
  gene_A = extract_first_match(`Alias(es) interactor A`, "(?<=uniprotkb:)[^|]+(?=\\(gene name\\))"),
  gene_B = extract_first_match(`Alias(es) interactor B`, "(?<=uniprotkb:)[^|]+(?=\\(gene name\\))"),
  taxid_A = extract_first_match(`Taxid interactor A`, "(?<=taxid:)\\d+"),
  taxid_B = extract_first_match(`Taxid interactor B`, "(?<=taxid:)\\d+"),
  type_A = `Type(s) interactor A`,
  type_B = `Type(s) interactor B`,
  interaction_type = `Interaction type(s)`,
  confidence = `Confidence value(s)`
)]

edge_dt <- edge_dt[
  taxid_A == "9606" & taxid_B == "9606" &
    type_A == 'psi-mi:"MI:0326"(protein)' &
    type_B == 'psi-mi:"MI:0326"(protein)' &
    !is.na(gene_A) & !is.na(gene_B) & gene_A != gene_B &
    interaction_type %in% c(
      'psi-mi:"MI:0407"(direct interaction)',
      'psi-mi:"MI:0915"(physical association)'
    )
]

edge_dt[, score := as.numeric(sub(".*intact-miscore:([0-9.]+).*", "\\1", confidence))]
edge_dt[!grepl("intact-miscore:", confidence), score := NA_real_]
edge_dt[, protein1 := pmin(gene_A, gene_B)]
edge_dt[, protein2 := pmax(gene_A, gene_B)]
edge_dt <- unique(edge_dt[, .(protein1, protein2, score)], by = c("protein1", "protein2"))

setorder(edge_dt, protein1, protein2)

write.csv(as.data.frame(edge_dt), output_path, row.names = TRUE)
cat("Filtered human physical IntAct edges:", nrow(edge_dt), "\n")
cat("Written:", output_path, "\n")
