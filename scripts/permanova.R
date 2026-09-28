args <- commandArgs(trailingOnly = TRUE)

count_file <- args[1]
metadata_file <- args[2]
output_file <- args[3]

library(vegan)

counts <- read.table(
  count_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE
)

metadata <- read.table(
  metadata_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

rownames(counts) <- counts$bin_id
counts$bin_id <- NULL

counts_t <- t(counts)

# Match metadata order to abundance matrix
metadata <- metadata[
  match(rownames(counts_t), metadata$sample_id),
]

if (any(is.na(metadata$sample_id))) {
  stop("Some samples are missing from metadata.")
}

# Check group replication
group_counts <- table(metadata$group)

if (length(group_counts) < 2 || any(group_counts < 2)) {

  message_text <- paste0(
    "PERMANOVA not performed.\n",
    "Reason: each comparison group requires biological replication.\n",
    "Group counts: ",
    paste(
      names(group_counts),
      group_counts,
      sep = "=",
      collapse = ", "
    ),
    "\n"
  )

  writeLines(
    message_text,
    output_file
  )

} else {

  bray <- vegdist(
    counts_t,
    method = "bray"
  )

  result <- adonis2(
    bray ~ group,
    data = metadata,
    permutations = 999
  )

  capture.output(
    result,
    file = output_file
  )
}
