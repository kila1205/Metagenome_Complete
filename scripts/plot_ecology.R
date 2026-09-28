#!/usr/bin/env Rscript

# ============================================================
# Andaliman WGS - Step 8 Ecology Plotting
#
# Generic plotting script for:
#   1. Alpha diversity
#   2. PCoA / beta diversity
#   3. Taxonomic composition
#   4. Ecological network
#
# No sample names or project-specific paths are hard-coded.
# ============================================================

suppressPackageStartupMessages({
    library(ggplot2)
    library(igraph)
})

# ------------------------------------------------------------
# Parse command-line arguments
# ------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(name, default = NULL) {
    pos <- match(name, args)

    if (!is.na(pos) && pos < length(args)) {
        return(args[pos + 1])
    }

    return(default)
}

metadata_file <- get_arg("--metadata")
alpha_file    <- get_arg("--alpha")
pcoa_file     <- get_arg("--pcoa")
relative_file <- get_arg("--relative")
taxonomy_file <- get_arg("--taxonomy")
edges_file    <- get_arg("--edges")
nodes_file    <- get_arg("--nodes")

outdir        <- get_arg("--outdir", "plots")
group_column  <- get_arg("--group_column", "group")
tax_rank      <- get_arg("--tax_rank", "genus")
top_taxa      <- as.integer(get_arg("--top_taxa", "20"))


# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------

read_tsv <- function(file) {

    if (is.null(file) || !file.exists(file)) {
        return(NULL)
    }

    read.delim(
        file,
        header = TRUE,
        sep = "\t",
        check.names = FALSE,
        stringsAsFactors = FALSE,
        comment.char = "",
        quote = ""
    )
}


find_column <- function(df, candidates) {

    if (is.null(df)) {
        return(NULL)
    }

    lower_names <- tolower(colnames(df))

    for (candidate in candidates) {

        hit <- which(lower_names == tolower(candidate))

        if (length(hit) > 0) {
            return(colnames(df)[hit[1]])
        }
    }

    return(NULL)
}


dir.create(
    outdir,
    recursive = TRUE,
    showWarnings = FALSE
)


# ============================================================
# 1. METADATA
# ============================================================

metadata <- read_tsv(metadata_file)

sample_column <- NULL

if (!is.null(metadata)) {

    sample_column <- find_column(
        metadata,
        c(
            "sample_id",
            "sample",
            "sampleid",
            "id"
        )
    )

    if (is.null(sample_column)) {
        sample_column <- colnames(metadata)[1]
    }

    message(
        "Metadata sample column: ",
        sample_column
    )

    if (group_column %in% colnames(metadata)) {

        message(
            "Grouping samples using metadata column: ",
            group_column
        )

    } else {

        message(
            "Metadata column '",
            group_column,
            "' not found. Plots will use sample identity only."
        )
    }
}


# ============================================================
# 2. ALPHA DIVERSITY
# ============================================================

alpha <- read_tsv(alpha_file)

if (!is.null(alpha)) {

    message("Creating alpha-diversity plots...")

    alpha_sample_column <- find_column(
        alpha,
        c(
            "sample_id",
            "sample",
            "sampleid",
            "id"
        )
    )

    if (is.null(alpha_sample_column)) {
        alpha_sample_column <- colnames(alpha)[1]
    }


    # Add metadata group if possible
    if (
        !is.null(metadata) &&
        group_column %in% colnames(metadata)
    ) {

        metadata_tmp <- metadata[
            ,
            c(sample_column, group_column),
            drop = FALSE
        ]

        colnames(metadata_tmp)[1] <- alpha_sample_column

        alpha <- merge(
            alpha,
            metadata_tmp,
            by = alpha_sample_column,
            all.x = TRUE,
            sort = FALSE
        )
    }


    excluded_columns <- c(
        alpha_sample_column,
        group_column
    )

    metric_columns <- setdiff(
        colnames(alpha),
        excluded_columns
    )

    metric_columns <- metric_columns[
        sapply(
            alpha[metric_columns],
            is.numeric
        )
    ]


    for (metric in metric_columns) {

        output_file <- file.path(
            outdir,
            paste0(
                "alpha_",
                tolower(gsub("[^A-Za-z0-9]+", "_", metric)),
                ".pdf"
            )
        )


        pdf(
            output_file,
            width = 8,
            height = 6
        )


        if (
            group_column %in% colnames(alpha) &&
            length(unique(na.omit(alpha[[group_column]]))) > 1
        ) {

            p <- ggplot(
                alpha,
                aes_string(
                    x = group_column,
                    y = metric
                )
            ) +
                geom_boxplot(
                    outlier.shape = NA
                ) +
                geom_jitter(
                    width = 0.15,
                    size = 3
                ) +
                labs(
                    title = paste(
                        metric,
                        "diversity"
                    ),
                    x = group_column,
                    y = metric
                ) +
                theme_bw(base_size = 12)

        } else {

            p <- ggplot(
                alpha,
                aes_string(
                    x = alpha_sample_column,
                    y = metric
                )
            ) +
                geom_col() +
                labs(
                    title = paste(
                        metric,
                        "diversity"
                    ),
                    x = "Sample",
                    y = metric
                ) +
                theme_bw(base_size = 12) +
                theme(
                    axis.text.x = element_text(
                        angle = 45,
                        hjust = 1
                    )
                )
        }


        print(p)

        dev.off()
    }
}


# ============================================================
# 3. PCoA
# ============================================================

pcoa <- read_tsv(pcoa_file)

if (!is.null(pcoa)) {

    message("Creating PCoA plot...")

    pcoa_sample_column <- find_column(
        pcoa,
        c(
            "sample_id",
            "sample",
            "sampleid",
            "id"
        )
    )

    if (is.null(pcoa_sample_column)) {
        pcoa_sample_column <- colnames(pcoa)[1]
    }


    numeric_columns <- colnames(pcoa)[
        sapply(
            pcoa,
            is.numeric
        )
    ]


    if (length(numeric_columns) >= 2) {

        axis1 <- numeric_columns[1]
        axis2 <- numeric_columns[2]


        if (
            !is.null(metadata) &&
            group_column %in% colnames(metadata)
        ) {

            metadata_tmp <- metadata[
                ,
                c(sample_column, group_column),
                drop = FALSE
            ]

            colnames(metadata_tmp)[1] <- pcoa_sample_column

            pcoa <- merge(
                pcoa,
                metadata_tmp,
                by = pcoa_sample_column,
                all.x = TRUE,
                sort = FALSE
            )
        }


        if (group_column %in% colnames(pcoa)) {

            p <- ggplot(
                pcoa,
                aes_string(
                    x = axis1,
                    y = axis2,
                    color = group_column
                )
            ) +
                geom_point(size = 4) +
                labs(
                    title = "PCoA - Bray-Curtis",
                    color = group_column
                )

        } else {

            p <- ggplot(
                pcoa,
                aes_string(
                    x = axis1,
                    y = axis2
                )
            ) +
                geom_point(size = 4)
        }


        p <- p +
            geom_text(
                aes_string(
                    label = pcoa_sample_column
                ),
                nudge_y = 0.02,
                check_overlap = TRUE
            ) +
            theme_bw(base_size = 12)


        ggsave(
            filename = file.path(
                outdir,
                "pcoa_bray_curtis.pdf"
            ),
            plot = p,
            width = 8,
            height = 6
        )
    }
}


# ============================================================
# 4. TAXONOMIC COMPOSITION
# ============================================================

relative <- read_tsv(relative_file)
taxonomy <- read_tsv(taxonomy_file)

if (
    !is.null(relative) &&
    !is.null(taxonomy) &&
    !is.null(metadata)
) {

    message(
        "Creating taxonomic composition plot..."
    )


    relative_id <- colnames(relative)[1]
    taxonomy_id <- colnames(taxonomy)[1]


    rank_column <- find_column(
        taxonomy,
        c(tax_rank)
    )


    if (!is.null(rank_column)) {

        taxonomy_tmp <- taxonomy[
            ,
            c(
                taxonomy_id,
                rank_column
            ),
            drop = FALSE
        ]

        colnames(taxonomy_tmp)[1] <- relative_id


        merged <- merge(
            relative,
            taxonomy_tmp,
            by = relative_id,
            all.x = TRUE
        )


        merged[[rank_column]][
            is.na(merged[[rank_column]]) |
            merged[[rank_column]] == ""
        ] <- "Unclassified"


        sample_names <- intersect(
            colnames(relative),
            metadata[[sample_column]]
        )


        if (length(sample_names) > 0) {

            long_list <- list()

            counter <- 1


            for (sample_name in sample_names) {

                long_list[[counter]] <- data.frame(
                    sample_id = sample_name,
                    taxon = merged[[rank_column]],
                    abundance = as.numeric(
                        merged[[sample_name]]
                    ),
                    stringsAsFactors = FALSE
                )

                counter <- counter + 1
            }


            long <- do.call(
                rbind,
                long_list
            )


            long$abundance[
                is.na(long$abundance)
            ] <- 0


            totals <- aggregate(
                abundance ~ taxon,
                data = long,
                FUN = sum
            )


            totals <- totals[
                order(
                    totals$abundance,
                    decreasing = TRUE
                ),
            ]


            keep_taxa <- head(
                totals$taxon,
                top_taxa
            )


            long$taxon_plot <- ifelse(
                long$taxon %in% keep_taxa,
                long$taxon,
                "Other"
            )


            collapsed <- aggregate(
                abundance ~ sample_id + taxon_plot,
                data = long,
                FUN = sum
            )


            p <- ggplot(
                collapsed,
                aes(
                    x = sample_id,
                    y = abundance,
                    fill = taxon_plot
                )
            ) +
                geom_col(
                    position = "stack"
                ) +
                labs(
                    title = paste(
                        "Taxonomic composition -",
                        rank_column
                    ),
                    x = "Sample",
                    y = "Relative abundance",
                    fill = rank_column
                ) +
                theme_bw(base_size = 12) +
                theme(
                    axis.text.x = element_text(
                        angle = 45,
                        hjust = 1
                    )
                )


            ggsave(
                filename = file.path(
                    outdir,
                    paste0(
                        "taxonomic_composition_",
                        tolower(rank_column),
                        ".pdf"
                    )
                ),
                plot = p,
                width = 10,
                height = 7
            )
        }
    } else {

        message(
            "Taxonomic rank '",
            tax_rank,
            "' was not found. Composition plot skipped."
        )
    }
}


# ============================================================
# 5. ECOLOGICAL NETWORK
# ============================================================

edges <- read_tsv(edges_file)
nodes <- read_tsv(nodes_file)

if (
    !is.null(edges) &&
    nrow(edges) > 0 &&
    ncol(edges) >= 2
) {

    message(
        "Creating ecological network plot..."
    )


    source_column <- find_column(
        edges,
        c(
            "source",
            "from",
            "node1",
            "taxon1"
        )
    )

    target_column <- find_column(
        edges,
        c(
            "target",
            "to",
            "node2",
            "taxon2"
        )
    )


    if (is.null(source_column)) {
        source_column <- colnames(edges)[1]
    }

    if (is.null(target_column)) {
        target_column <- colnames(edges)[2]
    }


    edge_data <- edges[
        ,
        c(
            source_column,
            target_column
        ),
        drop = FALSE
    ]


    graph <- graph_from_data_frame(
        edge_data,
        directed = FALSE
    )


    pdf(
        file.path(
            outdir,
            "ecological_network.pdf"
        ),
        width = 10,
        height = 10
    )


    plot(
        graph,
        vertex.size = 8,
        vertex.label.cex = 0.6,
        edge.width = 1,
        main = "Ecological association network"
    )


    dev.off()

} else {

    message(
        "No network edges detected. Network plot skipped."
    )
}


message(
    "Ecology plots completed successfully."
)
