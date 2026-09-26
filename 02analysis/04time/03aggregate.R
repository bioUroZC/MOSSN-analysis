library(data.table)

BASE_DIR <- "/proj/c.zihao/work1/02analysis/04time"
RESULTS_DIR <- file.path(BASE_DIR, "results")

read_timing_files <- function(pattern, required = TRUE) {
  files <- sort(list.files(RESULTS_DIR, pattern = pattern, full.names = TRUE))
  if (length(files) == 0L) {
    if (required) stop("No result files matched: ", pattern)
    return(data.table())
  }
  rbindlist(lapply(files, fread), fill = TRUE)
}

expected_grid <- data.table(
  n_samples = c(10L, 20L, 50L, 100L, 10L, 10L, 10L),
  n_edges = c(10000L, 10000L, 10000L, 10000L,
              20000L, 50000L, 100000L)
)

validate_timing_table <- function(timing_table, method_name) {
  required_columns <- c("method", "n_samples", "n_edges", "rep")
  if (!all(required_columns %in% names(timing_table))) {
    stop("Timing files for ", method_name, " lack required columns.")
  }
  if (!all(timing_table$method == method_name)) {
    stop("Timing files for ", method_name, " contain an unexpected method label.")
  }

  observed_grid <- unique(timing_table[, .(n_samples, n_edges)])
  if (nrow(fsetdiff(expected_grid, observed_grid)) > 0L ||
      nrow(fsetdiff(observed_grid, expected_grid)) > 0L) {
    stop("Timing files for ", method_name,
         " do not contain exactly the seven expected benchmark grids.")
  }

  repeats <- timing_table[, .N, by = .(n_samples, n_edges)]
  if (any(repeats$N != 3L)) {
    stop("Timing files for ", method_name,
         " must contain exactly three repeats per benchmark grid.")
  }
}

method_patterns <- c(
  SWEET = "^timing_SWEET_S[0-9]{3}_E[0-9]{6}\\.csv$",
  Patkar = "^timing_Patkar_S[0-9]{3}_E[0-9]{6}\\.csv$",
  PPIXpress = "^timing_PPIXpress_S[0-9]{3}_E[0-9]{6}\\.csv$",
  Proteinarium = "^timing_Proteinarium_S[0-9]{3}_E[0-9]{6}\\.csv$",
  MOSSN = "^timing_MOSSN_S[0-9]{3}_E[0-9]{6}\\.csv$",
  SSN = "^timing_pairedSSN_S[0-9]{3}_E[0-9]{6}\\.csv$"
)

py_tables <- vector("list", length(method_patterns))
names(py_tables) <- names(method_patterns)
for (method_name in names(method_patterns)) {
  method_table <- read_timing_files(method_patterns[[method_name]])
  validate_timing_table(method_table, method_name)
  py_tables[[method_name]] <- method_table
}
py <- rbindlist(py_tables, fill = TRUE)
setorderv(py, c("n_samples", "n_edges", "method", "rep"))

lioness <- read_timing_files("^timing_LIONESS_S[0-9]{3}_E[0-9]{6}\\.csv$")
validate_timing_table(lioness, "LIONESS")

fwrite(py, file.path(RESULTS_DIR, "timing_results_py.csv"))
fwrite(lioness, file.path(RESULTS_DIR, "timing_results_LIONESS.csv"))

cat("Saved:\n")
cat("  ", file.path(RESULTS_DIR, "timing_results_py.csv"), "\n", sep = "")
cat("  ", file.path(RESULTS_DIR, "timing_results_LIONESS.csv"), "\n", sep = "")
