args <- commandArgs(trailingOnly = TRUE)

relative_file <- args[1]
taxonomy_file <- args[2]
edges_file <- args[3]
nodes_file <- args[4]
message_file <- args[5]

relative <- read.table(
  relative_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

taxonomy <- read.table(
  taxonomy_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

rownames(relative) <- relative$bin_id
relative$bin_id <- NULL

n_samples <- ncol(relative)

if (n_samples < 3) {

  writeLines(
    paste0(
      "Ecological network not performed.\n",
      "Reason: correlation-based network analysis requires at least 3 samples.\n",
      "Current number of samples: ",
      n_samples,
      "\n"
    ),
    message_file
  )

  write.table(
    data.frame(
      source = character(),
      target = character(),
      correlation = numeric()
    ),
    edges_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  write.table(
    taxonomy,
    nodes_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

} else {

  # Correlation among bins/taxa across samples
  cor_matrix <- cor(
    t(relative),
    method = "spearman",
    use = "pairwise.complete.obs"
  )

  threshold <- 0.6

  edges <- data.frame(
    source = character(),
    target = character(),
    correlation = numeric()
  )

  taxa <- rownames(cor_matrix)

  for (i in 1:(nrow(cor_matrix) - 1)) {
    for (j in (i + 1):ncol(cor_matrix)) {

      r <- cor_matrix[i, j]

      if (!is.na(r) && abs(r) >= threshold) {

        edges <- rbind(
          edges,
          data.frame(
            source = taxa[i],
            target = taxa[j],
            correlation = r
          )
        )
      }
    }
  }

  write.table(
    edges,
    edges_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  write.table(
    taxonomy,
    nodes_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  writeLines(
    paste0(
      "Ecological network completed.\n",
      "Samples: ",
      n_samples,
      "\n",
      "Edges retained: ",
      nrow(edges),
      "\n",
      "Spearman threshold: |rho| >= ",
      threshold,
      "\n"
    ),
    message_file
  )
}
