rm(list = ls())

library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)

setwd("/proj/c.zihao/work1/1NT/9analysis/")

input_file <- "atlas_all.csv"
min_cancers <- 7
cancers <- c("BLCA", "BRCA", "CRC", "ESCA", "HNSC",
             "KIRC", "LIHC", "LUAD", "LUSC", "PRAD", "STAD")

atlas_all <- read.csv(
  input_file,
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)

atlas_sig <- atlas_all %>%
  dplyr::filter(direction %in% c("gain", "loss")) %>%
  dplyr::mutate(
    cancer = factor(cancer, levels = cancers),
    direction = factor(direction, levels = c("gain", "loss"))
  )

if (nrow(atlas_sig) == 0) {
  stop("No gain/loss links found in atlas_all.csv")
}

#=============================================================
# Step 2: define recurrently gained/lost links across cancers
#=============================================================

recurrent_tbl <- atlas_sig %>%
  dplyr::group_by(link) %>%
  dplyr::summarise(
    n_gain = base::sum(direction == "gain"),
    n_loss = base::sum(direction == "loss"),
    n_cancers = dplyr::n_distinct(cancer),
    median_delta_across_cancers = stats::median(delta_median, na.rm = TRUE),
    min_delta = min(delta_median, na.rm = TRUE),
    max_delta = max(delta_median, na.rm = TRUE),
    cancers_gain = paste(sort(unique(as.character(cancer[direction == "gain"]))), collapse = ","),
    cancers_loss = paste(sort(unique(as.character(cancer[direction == "loss"]))), collapse = ","),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    dominant_direction = dplyr::case_when(
      n_gain > n_loss ~ "gain",
      n_loss > n_gain ~ "loss",
      TRUE ~ "mixed"
    ),
    recurrent_count = pmax(n_gain, n_loss),
    consistency = dplyr::case_when(
      n_gain > 0 & n_loss == 0 ~ "ALL_GAIN",
      n_loss > 0 & n_gain == 0 ~ "ALL_LOSS",
      n_gain >= min_cancers & n_gain > n_loss ~ "MAJORITY_GAIN",
      n_loss >= min_cancers & n_loss > n_gain ~ "MAJORITY_LOSS",
      TRUE ~ "MIXED"
    ),
    recurrent_class = dplyr::case_when(
      n_gain >= min_cancers & n_gain > n_loss ~ "recurrently_gained",
      n_loss >= min_cancers & n_loss > n_gain ~ "recurrently_lost",
      TRUE ~ "not_recurrent"
    )
  ) %>%
  dplyr::arrange(
    dplyr::desc(recurrent_count),
    dplyr::desc(abs(median_delta_across_cancers))
  )

recurrent_links <- recurrent_tbl %>%
  dplyr::filter(recurrent_class != "not_recurrent")

recurrent_gain <- recurrent_tbl %>%
  dplyr::filter(recurrent_class == "recurrently_gained")

recurrent_loss <- recurrent_tbl %>%
  dplyr::filter(recurrent_class == "recurrently_lost")

write.csv(recurrent_links, "universal_recurrent_links.csv", row.names = FALSE)

base::message("Recurrent links kept: ", nrow(recurrent_links))
base::message("Recurrently gained: ", nrow(recurrent_gain))
base::message("Recurrently lost: ", nrow(recurrent_loss))

#=============================================================
# Step 3a: recurrence distribution
#=============================================================

recurrence_plot_df <- recurrent_tbl %>%
  dplyr::filter(recurrent_class != "not_recurrent") %>%
  dplyr::mutate(
    recurrence_bin = factor(
      recurrent_count,
      levels = seq(min_cancers, length(cancers), by = 1)
    ),
    recurrent_class = factor(
      recurrent_class,
      levels = c("recurrently_gained", "recurrently_lost")
    )
  )

recurrence_count_df <- recurrence_plot_df %>%
  dplyr::count(recurrence_bin, recurrent_class, name = "n_links")

p1 <- ggplot(
  recurrence_count_df,
  aes(x = recurrence_bin, fill = recurrent_class)
) +
  geom_col(aes(y = n_links), position = "stack", color = NA) +
  geom_text(
    aes(y = n_links, label = n_links),
    position = position_stack(vjust = 0.5),
    size = 3.2,
    color = "black"
  ) +
  scale_fill_manual(
    values = c(
      recurrently_gained = "#B2182B",
      recurrently_lost = "#2166AC"
    )
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  ) +
  labs(
    x = "Number of cancer types showing recurrent gain/loss",
    y = "Number of links",
    fill = "Class"
  )

pdf("recurrence_distribution.pdf", height = 5, width = 8)
print(p1)
dev.off()
