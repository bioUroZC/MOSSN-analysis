input_file <- "/proj/c.zihao/work1/01NT/03GEO/02eLUAD/eLUAD_exprSet_filtered.csv"
output_file <- "/proj/c.zihao/work1/01NT/03GEO/02eLUAD/eLUAD_normal_exprSet_filtered.csv"

expression <- read.csv(
  input_file,
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

normal_columns <- grep("_Normal$", colnames(expression), value = TRUE)
if (length(normal_columns) == 0L) {
  stop("No adjacent-normal samples ending in '_Normal' were found in: ", input_file)
}

normal_expression <- expression[, normal_columns, drop = FALSE]
write.csv(normal_expression, output_file, quote = FALSE)

message(
  "Wrote ", nrow(normal_expression), " genes x ", ncol(normal_expression),
  " cohort-matched normal samples to ", output_file
)
