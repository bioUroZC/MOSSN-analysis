library(data.table)

zip_path <- "/proj/c.zihao/work1/1NT/2intact/intact.zip"
run_full <- FALSE

normalize_intact_names <- function(dt) {
  setnames(dt, old = names(dt), new = sub("^#", "", names(dt)))
  dt
}

read_intact_from_zip <- function(zip_path, member = "intact.txt", nrows = Inf) {
  if (!file.exists(zip_path)) {
    stop("Zip file not found: ", zip_path)
  }

  members <- unzip(zip_path, list = TRUE)$Name
  if (!member %in% members) {
    stop("File not found in zip: ", member)
  }

  dt <- fread(
    cmd = sprintf("unzip -p %s %s", shQuote(zip_path), shQuote(member)),
    sep = "\t",
    header = TRUE,
    quote = "",
    fill = TRUE,
    nrows = nrows,
    showProgress = TRUE
  )

  normalize_intact_names(dt)
}

preview_intact_from_zip <- function(zip_path, member = "intact.txt", nrows = 5) {
  con <- pipe(
    sprintf(
      "unzip -p %s %s | head -n %d",
      shQuote(zip_path),
      shQuote(member),
      as.integer(nrows) + 1L
    ),
    open = "rt"
  )
  on.exit(close(con), add = TRUE)

  dt <- read.delim(
    con,
    sep = "\t",
    header = TRUE,
    check.names = FALSE
  )

  normalize_intact_names(as.data.table(dt))
}

extract_first_match <- function(x, pattern) {
  hit <- regexpr(pattern, x, perl = TRUE)
  out <- rep(NA_character_, length(x))
  ok <- hit > 0
  out[ok] <- regmatches(x, hit)
  out
}

preview_dt <- preview_intact_from_zip(zip_path, member = "intact.txt", nrows = 5)
preview_cols <- intersect(
  c(
    "ID(s) interactor A",
    "ID(s) interactor B",
    "Alias(es) interactor A",
    "Alias(es) interactor B",
    "Taxid interactor A",
    "Taxid interactor B",
    "Interaction type(s)"
  ),
  names(preview_dt)
)
print(preview_dt[, ..preview_cols])

cat("Preview rows:", nrow(preview_dt), "\n")






intact_dt <- read_intact_from_zip(zip_path, member = "intact.txt")

head(intact_dt)




edge_dt <- intact_dt[, .(
    protein_A = extract_first_match(`ID(s) interactor A`, "(?<=uniprotkb:)[^|]+"),
    protein_B = extract_first_match(`ID(s) interactor B`, "(?<=uniprotkb:)[^|]+"),
    gene_A = extract_first_match(`Alias(es) interactor A`, "(?<=uniprotkb:)[^|]+(?=\\(gene name\\))"),
    gene_B = extract_first_match(`Alias(es) interactor B`, "(?<=uniprotkb:)[^|]+(?=\\(gene name\\))"),
    taxid_A = extract_first_match(`Taxid interactor A`, "(?<=taxid:)\\d+"),
    taxid_B = extract_first_match(`Taxid interactor B`, "(?<=taxid:)\\d+"),
    type_A = `Type(s) interactor A`,
    type_B = `Type(s) interactor B`,
    interaction_type = `Interaction type(s)`,
    detection_method = `Interaction detection method(s)`,
    source_db = `Source database(s)`,
    interaction_id = `Interaction identifier(s)`,
    confidence = `Confidence value(s)`
  )]

head(edge_dt)


table(edge_dt$type_A)

edge_dt <- edge_dt[
  taxid_A == "9606" & taxid_B == "9606" &
    type_A == 'psi-mi:"MI:0326"(protein)' &
    type_B == 'psi-mi:"MI:0326"(protein)'
]



edge_dt <- edge_dt[!is.na(gene_A) & !is.na(gene_B)]

edge_dt <- edge_dt[gene_A != gene_B]

table(edge_dt$interaction_type)


edge_dt <- edge_dt[
  interaction_type %in% c(
    'psi-mi:"MI:0407"(direct interaction)',
    'psi-mi:"MI:0915"(physical association)'
  )
]


edge_dt[, miscore := as.numeric(sub(
  ".*intact-miscore:([0-9.]+).*",
  "\\1",
  confidence
))]
edge_dt[!grepl("intact-miscore:", confidence), miscore := NA_real_]


edge_dt[, gene_min := pmin(gene_A, gene_B)]
edge_dt[, gene_max := pmax(gene_A, gene_B)]

edge_dt <- unique(edge_dt, by = c("gene_min", "gene_max"))

names(edge_dt)

edge_dt <- subset(edge_dt, select = c("miscore", "gene_min", "gene_max"))

names(edge_dt) <- c("score", "protein1", "protein2")

cat("Filtered human PPIs:", nrow(edge_dt), "\n")
write.csv(edge_dt, file = 'intact_link.csv')

