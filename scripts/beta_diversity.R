args <- commandArgs(trailingOnly = TRUE)

count_file <- args[1]
distance_file <- args[2]
coordinate_file <- args[3]

library(vegan)

counts <- read.table(
  count_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE
)

rownames(counts) <- counts$bin_id
counts$bin_id <- NULL

counts_t <- t(counts)

# Bray-Curtis distance
bray <- vegdist(
  counts_t,
  method = "bray"
)

# Save distance matrix
bray_matrix <- as.matrix(bray)

write.table(
  bray_matrix,
  distance_file,
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

# Number of samples
n_samples <- nrow(counts_t)

if (n_samples < 2) {
  stop("Beta diversity requires at least 2 samples.")
}

# PCoA can have at most n_samples - 1 axes
k_axes <- min(2, n_samples - 1)

pcoa <- cmdscale(
  bray,
  k = k_axes,
  eig = TRUE
)

if (k_axes == 1) {

  coordinates <- data.frame(
    sample_id = rownames(pcoa$points),
    PCoA1 = pcoa$points[, 1],
    PCoA2 = 0
  )

} else {

  coordinates <- data.frame(
    sample_id = rownames(pcoa$points),
    PCoA1 = pcoa$points[, 1],
    PCoA2 = pcoa$points[, 2]
  )
}

write.table(
  coordinates,
  coordinate_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
