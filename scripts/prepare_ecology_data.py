#!/usr/bin/env python3

import pandas as pd
import sys
import os

abundance_file = sys.argv[1]
taxonomy_file = sys.argv[2]
metadata_file = sys.argv[3]
outdir = sys.argv[4]

os.makedirs(outdir, exist_ok=True)

# Read input tables
abundance = pd.read_csv(abundance_file, sep="\t")
taxonomy = pd.read_csv(taxonomy_file, sep="\t")
metadata = pd.read_csv(metadata_file, sep="\t")

# Standardize feature/bin identifier columns.
# This allows upstream tables to use names such as:
# bin, bin_id, feature_id, etc.
abundance = abundance.rename(
    columns={abundance.columns[0]: "bin_id"}
)

taxonomy = taxonomy.rename(
    columns={taxonomy.columns[0]: "bin_id"}
)

# Remove unmapped row if present
abundance = abundance[
    abundance["bin_id"].astype(str) != "unmapped"
].copy()

# Keep only bins found in taxonomy table
merged = abundance.merge(
    taxonomy,
    on="bin_id",
    how="inner"
)

if merged.empty:
    raise ValueError(
        "No matching bin IDs were found between abundance and taxonomy tables."
    )

# Taxonomic columns
taxonomy_cols = {
    "bin_id",
    "domain",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species"
}

# Sample columns come from abundance table
sample_cols = [
    col for col in abundance.columns
    if col not in taxonomy_cols and col != "bin_id"
]

if not sample_cols:
    raise ValueError(
        "No sample abundance columns were detected."
    )

# Check metadata
if "sample_id" not in metadata.columns:
    raise ValueError(
        "Metadata must contain a 'sample_id' column."
    )

metadata_samples = set(
    metadata["sample_id"].astype(str)
)

abundance_samples = set(
    map(str, sample_cols)
)

missing_metadata = abundance_samples - metadata_samples

if missing_metadata:
    raise ValueError(
        f"Samples missing from metadata: {sorted(missing_metadata)}"
    )

# Count matrix
count_matrix = merged[
    ["bin_id"] + sample_cols
].copy()

# Relative abundance
relative = count_matrix.copy()

for sample in sample_cols:

    total = relative[sample].sum()

    if total > 0:
        relative[sample] = relative[sample] / total
    else:
        relative[sample] = 0

# Make sure expected taxonomy columns exist.
# Missing classifications are allowed.
for col in [
    "domain",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species"
]:
    if col not in merged.columns:
        merged[col] = ""

# Taxonomy table
taxonomy_clean = merged[
    [
        "bin_id",
        "domain",
        "phylum",
        "class",
        "order",
        "family",
        "genus",
        "species"
    ]
].copy()

# Save outputs
count_matrix.to_csv(
    f"{outdir}/count_matrix.tsv",
    sep="\t",
    index=False
)

relative.to_csv(
    f"{outdir}/relative_abundance.tsv",
    sep="\t",
    index=False
)

taxonomy_clean.to_csv(
    f"{outdir}/taxonomy_clean.tsv",
    sep="\t",
    index=False
)

metadata.to_csv(
    f"{outdir}/metadata_clean.tsv",
    sep="\t",
    index=False
)

print(f"Prepared {len(count_matrix)} bins")
print(f"Samples: {', '.join(sample_cols)}")
