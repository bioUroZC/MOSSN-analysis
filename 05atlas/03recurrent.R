rm(list = ls())

base_dir <- "/proj/c.zihao/work1/05atlas"
min_cancers <- 7

atlas_all <- read.csv(file.path(base_dir, "atlas_all.csv"), row.names = 1)
atlas_sig <- atlas_all[atlas_all$direction %in% c("gained", "lost"), ]

links <- sort(unique(atlas_sig$link), method = "radix")
rows_by_link <- split(seq_len(nrow(atlas_sig)), factor(atlas_sig$link, levels = links))

n_link <- length(links)
n_gained <- numeric(n_link)
n_lost <- numeric(n_link)
n_cancers <- numeric(n_link)
median_delta <- numeric(n_link)
min_delta <- numeric(n_link)
max_delta <- numeric(n_link)
cancers_gained <- character(n_link)
cancers_lost <- character(n_link)

for (i in seq_len(n_link)) {
  one <- atlas_sig[rows_by_link[[i]], ]
  is_gained <- one$direction == "gained"
  is_lost <- one$direction == "lost"

  n_gained[i] <- sum(is_gained)
  n_lost[i] <- sum(is_lost)
  n_cancers[i] <- length(unique(one$cancer))
  median_delta[i] <- median(one$delta_median)
  min_delta[i] <- min(one$delta_median)
  max_delta[i] <- max(one$delta_median)
  cancers_gained[i] <- paste(sort(unique(one$cancer[is_gained])), collapse = ",")
  cancers_lost[i] <- paste(sort(unique(one$cancer[is_lost])), collapse = ",")
}

recurrent_tbl <- data.frame(
  link = links,
  n_gained = n_gained,
  n_lost = n_lost,
  n_cancers = n_cancers,
  median_delta_across_cancers = median_delta,
  min_delta = min_delta,
  max_delta = max_delta,
  cancers_gained = cancers_gained,
  cancers_lost = cancers_lost
)

recurrent_tbl$dominant_direction <- "mixed"
recurrent_tbl$dominant_direction[n_gained > n_lost] <- "gained"
recurrent_tbl$dominant_direction[n_lost > n_gained] <- "lost"

recurrent_tbl$recurrent_count <- pmax(n_gained, n_lost)

recurrent_tbl$consistency <- "MIXED"
recurrent_tbl$consistency[n_lost >= min_cancers & n_lost > n_gained] <- "MAJORITY_LOST"
recurrent_tbl$consistency[n_gained >= min_cancers & n_gained > n_lost] <- "MAJORITY_GAINED"
recurrent_tbl$consistency[n_lost > 0 & n_gained == 0] <- "ALL_LOST"
recurrent_tbl$consistency[n_gained > 0 & n_lost == 0] <- "ALL_GAINED"

recurrent_tbl$recurrent_class <- "not_recurrent"
recurrent_tbl$recurrent_class[n_gained >= min_cancers & n_gained > n_lost] <- "recurrently_gained"
recurrent_tbl$recurrent_class[n_lost >= min_cancers & n_lost > n_gained] <- "recurrently_lost"

recurrent_tbl <- recurrent_tbl[order(-recurrent_tbl$recurrent_count,
                                     -abs(recurrent_tbl$median_delta_across_cancers)), ]

recurrent_links <- recurrent_tbl[recurrent_tbl$recurrent_class != "not_recurrent", ]
write.csv(recurrent_links, file.path(base_dir, "universal_recurrent_links.csv"),
          row.names = FALSE)

message("Recurrent links kept: ", nrow(recurrent_links))
message("Recurrently gained: ", sum(recurrent_links$recurrent_class == "recurrently_gained"))
message("Recurrently lost: ", sum(recurrent_links$recurrent_class == "recurrently_lost"))
