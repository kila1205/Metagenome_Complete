#!/usr/bin/env python3

import csv
import sys

input_file = sys.argv[1]
output_file = sys.argv[2]

ranks = [
    "domain",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species"
]

with open(input_file, newline="") as infile, \
     open(output_file, "w", newline="") as outfile:

    reader = csv.DictReader(infile, delimiter="\t")

    writer = csv.writer(outfile, delimiter="\t")
    writer.writerow(["bin_id"] + ranks)

    for row in reader:
        bin_id = row["user_genome"]
        classification = row["classification"]

        parts = classification.split(";")

        taxonomy = {}

        for part in parts:
            if "__" in part:
                prefix, name = part.split("__", 1)

                rank_map = {
                    "d": "domain",
                    "p": "phylum",
                    "c": "class",
                    "o": "order",
                    "f": "family",
                    "g": "genus",
                    "s": "species"
                }

                if prefix in rank_map:
                    taxonomy[rank_map[prefix]] = name

        writer.writerow([
            bin_id,
            taxonomy.get("domain", "Unclassified"),
            taxonomy.get("phylum", "Unclassified"),
            taxonomy.get("class", "Unclassified"),
            taxonomy.get("order", "Unclassified"),
            taxonomy.get("family", "Unclassified"),
            taxonomy.get("genus", "Unclassified"),
            taxonomy.get("species", "Unclassified")
        ])
