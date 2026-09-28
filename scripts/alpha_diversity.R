args <- commandArgs(trailingOnly = TRUE)

count_file <- args[1]
output_file <- args[2]

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

shannon <- diversity(counts_t, index = "shannon")
simpson <- diversity(counts_t, index = "simpson")

observed <- specnumber(counts_t)

chao <- estimateR(counts_t)
chao1 <- chao["S.chao1", ]

evenness <- shannon / log(observed)

result <- data.frame(
  sample_id = rownames(counts_t),
  observed_richness = observed,
  shannon = shannon,
  simpson = simpson,
  chao1 = chao1,
  evenness = evenness
)

write.table(
  result,
  output_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
