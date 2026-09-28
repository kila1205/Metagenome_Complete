#!/usr/bin/env python3

import argparse
import csv
import os
import re
from collections import defaultdict

import matplotlib.pyplot as plt

from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import cm
from reportlab.platypus import (
    SimpleDocTemplate,
    Paragraph,
    Spacer,
    Table,
    TableStyle,
    PageBreak,
    Image,
    KeepTogether,
)


# ============================================================
# Arguments
# ============================================================

parser = argparse.ArgumentParser()

parser.add_argument("--sample-id", required=True)
parser.add_argument("--raw-nanostats", required=True)
parser.add_argument("--filtered-nanostats", required=True)
parser.add_argument("--assembly-stats", required=True)
parser.add_argument("--quast-report", required=True)
parser.add_argument("--protein-stats", required=True)
parser.add_argument("--bin-stats", required=True)
parser.add_argument("--coverm", required=True)
parser.add_argument("--kraken-report", required=True)
parser.add_argument("--bracken", required=True)
parser.add_argument("--taxonomy", required=True)
parser.add_argument("--catalog-stats", required=True)
parser.add_argument("--kofam", required=True)
parser.add_argument("--eggnog", required=True)

parser.add_argument("--output-pdf", required=True)
parser.add_argument("--output-summary", required=True)
parser.add_argument("--figures-dir", required=True)

args = parser.parse_args()

os.makedirs(args.figures_dir, exist_ok=True)


# ============================================================
# Generic helpers
# ============================================================

def safe_float(value, default=0.0):
    try:
        return float(str(value).replace(",", ""))
    except Exception:
        return default


def safe_int(value, default=0):
    try:
        return int(float(str(value).replace(",", "")))
    except Exception:
        return default


def read_tsv(path, comments=False):
    rows = []

    with open(path, "r", encoding="utf-8", errors="replace") as handle:

        if comments:
            lines = [
                line for line in handle
                if line.strip() and not line.startswith("#")
            ]
        else:
            lines = [
                line for line in handle
                if line.strip()
            ]

    if not lines:
        return []

    reader = csv.DictReader(lines, delimiter="\t")

    for row in reader:
        rows.append(row)

    return rows


def find_first(directory_or_file, basename):
    if os.path.isfile(directory_or_file):
        return directory_or_file

    for root, _, files in os.walk(directory_or_file):
        if basename in files:
            return os.path.join(root, basename)

    return None


def human_number(number):
    number = float(number)

    if number >= 1_000_000_000:
        return f"{number / 1_000_000_000:.2f} Gb"
    if number >= 1_000_000:
        return f"{number / 1_000_000:.2f} Mb"
    if number >= 1_000:
        return f"{number / 1_000:.2f} kb"

    return f"{number:.0f}"


# ============================================================
# NanoPlot parsing
# ============================================================

def parse_nanostats(path):

    result = {}

    patterns = {
        "mean_length": r"Mean read length:\s*([\d,\.]+)",
        "mean_quality": r"Mean read quality:\s*([\d,\.]+)",
        "median_length": r"Median read length:\s*([\d,\.]+)",
        "median_quality": r"Median read quality:\s*([\d,\.]+)",
        "number_reads": r"Number of reads:\s*([\d,]+)",
        "read_length_n50": r"Read length N50:\s*([\d,\.]+)",
        "total_bases": r"Total bases:\s*([\d,]+)",
    }

    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        text = handle.read()

    for key, pattern in patterns.items():

        match = re.search(pattern, text)

        if match:
            value = match.group(1).replace(",", "")

            if key in {
                "mean_quality",
                "median_quality",
                "mean_length",
                "median_length",
                "read_length_n50",
            }:
                result[key] = safe_float(value)
            else:
                result[key] = safe_int(value)

    return result


# ============================================================
# SeqKit stats parsing
# ============================================================

def parse_seqkit_stats(path):

    rows = read_tsv(path)

    if not rows:
        return {}

    return rows[0]


# ============================================================
# QUAST parsing
# ============================================================

def parse_quast(path):

    report_file = find_first(path, "report.tsv")

    if not report_file:
        return {}

    result = {}

    with open(
        report_file,
        "r",
        encoding="utf-8",
        errors="replace"
    ) as handle:

        reader = csv.reader(handle, delimiter="\t")

        for row in reader:

            if len(row) >= 2:
                result[row[0].strip()] = row[1].strip()

    return result


# ============================================================
# Bin statistics
# ============================================================

bin_rows = read_tsv(args.bin_stats)

bin_count = len(bin_rows)

bin_total_length = sum(
    safe_float(row.get("length", 0))
    for row in bin_rows
)


# ============================================================
# CoverM parsing
# ============================================================

coverm_rows = read_tsv(args.coverm)

coverm_abundance = {}

for row in coverm_rows:

    if not row:
        continue

    first_col = list(row.keys())[0]
    bin_id = row.get(first_col, "")

    abundance_value = 0.0

    for key, value in row.items():

        normalized_key = re.sub(r"[^a-z0-9]+", "_", key.lower()).strip("_")

        if "relative_abundance" in normalized_key:
            abundance_value = safe_float(value)
            break

    coverm_abundance[bin_id] = abundance_value


# ============================================================
# GTDB taxonomy
# ============================================================

taxonomy_rows = read_tsv(args.taxonomy)

taxonomy_by_bin = {}

for row in taxonomy_rows:

    bin_id = row.get("bin_id", "")

    if bin_id:
        taxonomy_by_bin[bin_id] = row


# ============================================================
# Aggregate GTDB taxonomy + CoverM
# ============================================================

species_abundance = defaultdict(float)
species_bins = defaultdict(int)

for bin_id, taxonomy in taxonomy_by_bin.items():

    species = taxonomy.get("species", "").strip()

    if not species:
        species = "Unclassified"

    abundance = coverm_abundance.get(bin_id, 0.0)

    species_abundance[species] += abundance
    species_bins[species] += 1


gtdb_top = sorted(
    species_abundance.items(),
    key=lambda x: x[1],
    reverse=True
)


# ============================================================
# Kraken2
# ============================================================

kraken_top = []

with open(
    args.kraken_report,
    "r",
    encoding="utf-8",
    errors="replace"
) as handle:

    for line in handle:

        parts = line.rstrip("\n").split("\t")

        if len(parts) < 6:
            continue

        percentage = safe_float(parts[0])
        rank = parts[3].strip()
        taxon = parts[5].strip()

        if rank == "S":
            kraken_top.append(
                (taxon, percentage)
            )

kraken_top.sort(
    key=lambda x: x[1],
    reverse=True
)


# ============================================================
# Bracken
# ============================================================

bracken_rows = read_tsv(args.bracken)

bracken_top = []

for row in bracken_rows:

    name = (
        row.get("name")
        or row.get("taxonomy_name")
        or row.get("Name")
        or "Unknown"
    )

    fraction = (
        row.get("fraction_total_reads")
        or row.get("fraction")
        or row.get("new_est_reads")
        or 0
    )

    fraction = safe_float(fraction)

    # Bracken fraction_total_reads is normally 0-1.
    if fraction <= 1:
        fraction *= 100

    bracken_top.append(
        (name, fraction)
    )

bracken_top.sort(
    key=lambda x: x[1],
    reverse=True
)


# ============================================================
# Functional annotation counts
# ============================================================

# KofamScan detail-tsv can contain many candidate profile matches per gene.
# Rows beginning with "*" are the threshold-passing assignments reported by
# KofamScan. Count those separately from all candidate matches.
kofam_candidate_count = 0
kofam_accepted_count = 0
kofam_accepted_kos = set()
kofam_accepted_definitions = defaultdict(int)

with open(args.kofam, "r", encoding="utf-8", errors="replace") as handle:
    for raw_line in handle:
        if not raw_line.strip() or raw_line.lstrip().startswith("#"):
            continue

        kofam_candidate_count += 1
        accepted = raw_line.startswith("*") or raw_line.lstrip().startswith("*")
        if not accepted:
            continue

        kofam_accepted_count += 1
        fields = raw_line.lstrip("* \t").rstrip("\n").split("\t")
        if len(fields) < 2:
            fields = raw_line.lstrip("* \t").split(None, 5)

        if len(fields) >= 2:
            kofam_accepted_kos.add(fields[1].strip())
        if len(fields) >= 6:
            kofam_accepted_definitions[fields[5].strip().strip('"')] += 1

# eggNOG-mapper output: count annotation rows and useful functional fields.
eggnog_count = 0
eggnog_cog_count = 0
eggnog_ko_count = 0
eggnog_pathway_count = 0
eggnog_cog_categories = defaultdict(int)

with open(args.eggnog, "r", encoding="utf-8", errors="replace") as handle:
    header = None
    for raw_line in handle:
        line = raw_line.rstrip("\n")
        if line.startswith("#query"):
            header = line.lstrip("#").split("\t")
            continue
        if not line or line.startswith("#"):
            continue

        eggnog_count += 1
        fields = line.split("\t")
        if not header or len(fields) < len(header):
            continue
        row = dict(zip(header, fields))

        cog = row.get("COG_category", "-").strip()
        ko = row.get("KEGG_ko", "-").strip()
        pathway = row.get("KEGG_Pathway", "-").strip()

        if cog not in {"", "-"}:
            eggnog_cog_count += 1
            for category in cog:
                if category.isalpha():
                    eggnog_cog_categories[category] += 1
        if ko not in {"", "-"}:
            eggnog_ko_count += 1
        if pathway not in {"", "-"}:
            eggnog_pathway_count += 1

# Keep legacy variable name for downstream compatibility, but make it mean
# biologically accepted KofamScan assignments rather than all candidate rows.
kofam_count = kofam_accepted_count

# ============================================================
# Core QC/statistics
# ============================================================

raw_qc = parse_nanostats(args.raw_nanostats)
filtered_qc = parse_nanostats(args.filtered_nanostats)

assembly_stats = parse_seqkit_stats(args.assembly_stats)
protein_stats = parse_seqkit_stats(args.protein_stats)
catalog_stats = parse_seqkit_stats(args.catalog_stats)

quast = parse_quast(args.quast_report)


# ============================================================

# ============================================================
# REPORT PRESENTATION LAYER
# ============================================================

from collections import Counter
from reportlab.lib.enums import TA_CENTER, TA_LEFT, TA_JUSTIFY


# ============================================================
# Additional functional summaries for visualization
# ============================================================

# ------------------------------------------------------------
# KofamScan: accepted biological functions
# ------------------------------------------------------------

kofam_function_counts = Counter()

with open(
    args.kofam,
    "r",
    encoding="utf-8",
    errors="replace"
) as handle:

    for raw_line in handle:

        if (
            not raw_line.strip()
            or raw_line.lstrip().startswith("#")
        ):
            continue

        # Only threshold-passing KofamScan assignments.
        if not (
            raw_line.startswith("*")
            or raw_line.lstrip().startswith("*")
        ):
            continue

        cleaned = raw_line.lstrip("* \t").rstrip("\n")

        fields = cleaned.split("\t")

        if len(fields) < 6:
            fields = cleaned.split(None, 5)

        if len(fields) >= 6:

            definition = fields[5].strip().strip('"')

            if definition:
                kofam_function_counts[definition] += 1


# ------------------------------------------------------------
# eggNOG: COG and KEGG pathway summaries
# ------------------------------------------------------------

eggnog_cog_plot_counts = Counter()
eggnog_pathway_counts = Counter()

with open(
    args.eggnog,
    "r",
    encoding="utf-8",
    errors="replace"
) as handle:

    header = None

    for raw_line in handle:

        line = raw_line.rstrip("\n")

        if line.startswith("#query"):
            header = line.lstrip("#").split("\t")
            continue

        if not line or line.startswith("#"):
            continue

        fields = line.split("\t")

        if not header:
            continue

        row = dict(zip(header, fields))

        cog = row.get(
            "COG_category",
            "-"
        ).strip()

        if cog not in {"", "-"}:

            for category in cog:

                if category.isalpha():
                    eggnog_cog_plot_counts[category] += 1

        pathways = row.get(
            "KEGG_Pathway",
            "-"
        ).strip()

        if pathways not in {"", "-"}:

            for pathway in pathways.split(","):

                pathway = pathway.strip()

                # Prefer koXXXXX and ignore duplicated mapXXXXX entries.
                if re.fullmatch(
                    r"ko\d{5}",
                    pathway
                ):
                    eggnog_pathway_counts[pathway] += 1


# ============================================================
# Biological labels
# ============================================================

COG_LABELS = {
    "A": "RNA processing and modification",
    "B": "Chromatin structure and dynamics",
    "C": "Energy production and conversion",
    "D": "Cell cycle and division",
    "E": "Amino acid metabolism and transport",
    "F": "Nucleotide metabolism and transport",
    "G": "Carbohydrate metabolism and transport",
    "H": "Coenzyme metabolism",
    "I": "Lipid metabolism",
    "J": "Translation and ribosome biogenesis",
    "K": "Transcription",
    "L": "Replication, recombination and repair",
    "M": "Cell wall and membrane biogenesis",
    "N": "Cell motility",
    "O": "Protein turnover and chaperones",
    "P": "Inorganic ion transport and metabolism",
    "Q": "Secondary metabolite biosynthesis",
    "S": "Function unknown",
    "T": "Signal transduction",
    "U": "Intracellular trafficking and secretion",
    "V": "Defense mechanisms",
    "W": "Extracellular structures",
    "Z": "Cytoskeleton",
}


KEGG_PATHWAY_LABELS = {
    "ko01100": "Metabolic pathways",
    "ko01110": "Biosynthesis of secondary metabolites",
    "ko01120": "Microbial metabolism in diverse environments",
    "ko01130": "Biosynthesis of antibiotics",
    "ko02010": "ABC transporters",
    "ko02020": "Two-component system",
    "ko01230": "Biosynthesis of amino acids",
    "ko01200": "Carbon metabolism",
    "ko02024": "Quorum sensing",
    "ko00230": "Purine metabolism",
    "ko00240": "Pyrimidine metabolism",
    "ko00620": "Pyruvate metabolism",
    "ko00010": "Glycolysis / gluconeogenesis",
    "ko00520": "Amino sugar and nucleotide sugar metabolism",
    "ko00500": "Starch and sucrose metabolism",
    "ko03070": "Bacterial secretion system",
    "ko00270": "Cysteine and methionine metabolism",
    "ko00190": "Oxidative phosphorylation",
    "ko00630": "Glyoxylate and dicarboxylate metabolism",
    "ko00260": "Glycine, serine and threonine metabolism",
}


# ============================================================
# Figures
# ============================================================

qc_figure = os.path.join(
    args.figures_dir,
    "qc_comparison.png"
)

raw_reads = raw_qc.get(
    "number_reads",
    0
)

filtered_reads = filtered_qc.get(
    "number_reads",
    0
)

raw_bases = raw_qc.get(
    "total_bases",
    0
)

filtered_bases = filtered_qc.get(
    "total_bases",
    0
)


# ------------------------------------------------------------
# QC figure
# ------------------------------------------------------------

fig, ax = plt.subplots(
    figsize=(7.2, 4.3)
)

labels = [
    "Reads\n(millions)",
    "Bases\n(Gb)"
]

raw_values = [
    raw_reads / 1_000_000,
    raw_bases / 1_000_000_000
]

filtered_values = [
    filtered_reads / 1_000_000,
    filtered_bases / 1_000_000_000
]

x = range(
    len(labels)
)

ax.bar(
    [i - 0.18 for i in x],
    raw_values,
    width=0.36,
    label="Raw"
)

ax.bar(
    [i + 0.18 for i in x],
    filtered_values,
    width=0.36,
    label="Filtered"
)

ax.set_xticks(
    list(x)
)

ax.set_xticklabels(
    labels
)

ax.set_ylabel(
    "Amount of sequencing data"
)

ax.set_title(
    "Sequencing data before and after filtering"
)

ax.legend()

fig.tight_layout()

fig.savefig(
    qc_figure,
    dpi=180,
    bbox_inches="tight"
)

plt.close(fig)


# ------------------------------------------------------------
# Read-based taxonomy
# ------------------------------------------------------------

taxonomy_figure = os.path.join(
    args.figures_dir,
    "read_taxonomy_top_species.png"
)

if bracken_top:

    top = bracken_top[:15]

    names = [
        name
        for name, _ in reversed(top)
    ]

    values = [
        value
        for _, value in reversed(top)
    ]

    fig, ax = plt.subplots(
        figsize=(8.3, 6.0)
    )

    ax.barh(
        names,
        values
    )

    ax.set_xlabel(
        "Estimated relative abundance (%)"
    )

    ax.set_title(
        "Most abundant species in the read-based community profile"
    )

    fig.tight_layout()

    fig.savefig(
        taxonomy_figure,
        dpi=180,
        bbox_inches="tight"
    )

    plt.close(fig)


# ------------------------------------------------------------
# Genome-resolved taxonomy
# ------------------------------------------------------------

gtdb_figure = os.path.join(
    args.figures_dir,
    "bin_taxonomy_abundance.png"
)

if gtdb_top:

    top = gtdb_top[:15]

    names = [
        name
        for name, _ in reversed(top)
    ]

    values = [
        value
        for _, value in reversed(top)
    ]

    fig, ax = plt.subplots(
        figsize=(8.3, 6.0)
    )

    ax.barh(
        names,
        values
    )

    ax.set_xlabel(
        "CoverM relative abundance (%)"
    )

    ax.set_title(
        "Most abundant reconstructed genome populations"
    )

    fig.tight_layout()

    fig.savefig(
        gtdb_figure,
        dpi=180,
        bbox_inches="tight"
    )

    plt.close(fig)


# ------------------------------------------------------------
# KofamScan biological-function figure
# ------------------------------------------------------------

kofam_figure = os.path.join(
    args.figures_dir,
    "kofamscan_top_functions.png"
)

if kofam_function_counts:

    top_kofam = kofam_function_counts.most_common(12)

    function_names = []
    function_values = []

    for name, count in reversed(top_kofam):

        clean_name = name

        if len(clean_name) > 62:
            clean_name = (
                clean_name[:59]
                + "..."
            )

        function_names.append(
            clean_name
        )

        function_values.append(
            count
        )

    fig, ax = plt.subplots(
        figsize=(9.2, 6.2)
    )

    ax.barh(
        function_names,
        function_values
    )

    ax.set_xlabel(
        "Accepted KofamScan assignments"
    )

    ax.set_title(
        "Most frequently predicted protein functions"
    )

    fig.tight_layout()

    fig.savefig(
        kofam_figure,
        dpi=180,
        bbox_inches="tight"
    )

    plt.close(fig)


# ------------------------------------------------------------
# eggNOG COG figure
# ------------------------------------------------------------

cog_figure = os.path.join(
    args.figures_dir,
    "eggnog_cog_categories.png"
)

if eggnog_cog_plot_counts:

    top_cogs = sorted(
        eggnog_cog_plot_counts.items(),
        key=lambda x: x[1],
        reverse=True
    )[:15]

    cog_names = []
    cog_values = []

    for code, count in reversed(top_cogs):

        label = COG_LABELS.get(
            code,
            f"COG category {code}"
        )

        cog_names.append(
            label
        )

        cog_values.append(
            count
        )

    fig, ax = plt.subplots(
        figsize=(9.2, 6.5)
    )

    ax.barh(
        cog_names,
        cog_values
    )

    ax.set_xlabel(
        "Annotated proteins"
    )

    ax.set_title(
        "Major functional categories in the protein catalog"
    )

    fig.tight_layout()

    fig.savefig(
        cog_figure,
        dpi=180,
        bbox_inches="tight"
    )

    plt.close(fig)


# ------------------------------------------------------------
# eggNOG KEGG pathway figure
# ------------------------------------------------------------

pathway_figure = os.path.join(
    args.figures_dir,
    "eggnog_kegg_pathways.png"
)

if eggnog_pathway_counts:

    top_pathways = eggnog_pathway_counts.most_common(15)

    pathway_names = []
    pathway_values = []

    for pathway, count in reversed(top_pathways):

        label = KEGG_PATHWAY_LABELS.get(
            pathway,
            pathway
        )

        pathway_names.append(
            label
        )

        pathway_values.append(
            count
        )

    fig, ax = plt.subplots(
        figsize=(9.2, 6.5)
    )

    ax.barh(
        pathway_names,
        pathway_values
    )

    ax.set_xlabel(
        "Proteins associated with pathway"
    )

    ax.set_title(
        "Most represented KEGG pathways"
    )

    fig.tight_layout()

    fig.savefig(
        pathway_figure,
        dpi=180,
        bbox_inches="tight"
    )

    plt.close(fig)


# ============================================================
# Summary TSV
# ============================================================

summary_rows = [
    (
        "sample_id",
        args.sample_id
    ),
    (
        "raw_reads",
        raw_reads
    ),
    (
        "raw_bases",
        raw_bases
    ),
    (
        "filtered_reads",
        filtered_reads
    ),
    (
        "filtered_bases",
        filtered_bases
    ),
    (
        "assembly_sequences",
        assembly_stats.get(
            "num_seqs",
            ""
        )
    ),
    (
        "assembly_length",
        assembly_stats.get(
            "sum_len",
            ""
        )
    ),
    (
        "assembly_n50",
        assembly_stats.get(
            "N50",
            ""
        )
    ),
    (
        "predicted_proteins",
        protein_stats.get(
            "num_seqs",
            ""
        )
    ),
    (
        "candidate_bins",
        bin_count
    ),
    (
        "representative_proteins",
        catalog_stats.get(
            "num_seqs",
            ""
        )
    ),
    (
        "kofam_accepted_assignments",
        kofam_accepted_count
    ),
    (
        "kofam_unique_kos",
        len(kofam_accepted_kos)
    ),
    (
        "kofam_candidate_matches",
        kofam_candidate_count
    ),
    (
        "eggnog_annotation_rows",
        eggnog_count
    ),
    (
        "eggnog_rows_with_cog",
        eggnog_cog_count
    ),
    (
        "eggnog_rows_with_kegg_ko",
        eggnog_ko_count
    ),
    (
        "eggnog_rows_with_kegg_pathway",
        eggnog_pathway_count
    ),
]


with open(
    args.output_summary,
    "w",
    encoding="utf-8",
    newline=""
) as handle:

    writer = csv.writer(
        handle,
        delimiter="\t"
    )

    writer.writerow([
        "metric",
        "value"
    ])

    writer.writerows(
        summary_rows
    )


# ============================================================
# PDF styles
# ============================================================

styles = getSampleStyleSheet()


styles.add(
    ParagraphStyle(
        name="ReportTitle",
        parent=styles["Title"],
        fontSize=22,
        leading=27,
        alignment=TA_CENTER,
        spaceAfter=12,
    )
)


styles.add(
    ParagraphStyle(
        name="ReportSubtitle",
        parent=styles["BodyText"],
        fontSize=11,
        leading=15,
        alignment=TA_CENTER,
        spaceAfter=16,
    )
)


styles.add(
    ParagraphStyle(
        name="SectionTitle",
        parent=styles["Heading1"],
        fontSize=15,
        leading=19,
        alignment=TA_LEFT,
        spaceBefore=10,
        spaceAfter=9,
        keepWithNext=True,
    )
)


styles.add(
    ParagraphStyle(
        name="SubsectionTitle",
        parent=styles["Heading2"],
        fontSize=11.5,
        leading=15,
        alignment=TA_LEFT,
        spaceBefore=8,
        spaceAfter=6,
        keepWithNext=True,
    )
)


styles.add(
    ParagraphStyle(
        name="BodyReport",
        parent=styles["BodyText"],
        fontSize=9.5,
        leading=14,
        alignment=TA_JUSTIFY,
        spaceAfter=8,
    )
)


styles.add(
    ParagraphStyle(
        name="Note",
        parent=styles["BodyText"],
        fontSize=8.7,
        leading=12.5,
        alignment=TA_JUSTIFY,
        leftIndent=10,
        rightIndent=10,
        spaceBefore=4,
        spaceAfter=8,
    )
)


styles.add(
    ParagraphStyle(
        name="MetricNumber",
        parent=styles["BodyText"],
        fontSize=13,
        leading=16,
        alignment=TA_CENTER,
    )
)


styles.add(
    ParagraphStyle(
        name="MetricLabel",
        parent=styles["BodyText"],
        fontSize=8,
        leading=10,
        alignment=TA_CENTER,
    )
)


styles.add(
    ParagraphStyle(
        name="TableText",
        parent=styles["BodyText"],
        fontSize=8,
        leading=10,
        alignment=TA_LEFT,
    )
)


styles.add(
    ParagraphStyle(
        name="TableHeader",
        parent=styles["BodyText"],
        fontSize=8,
        leading=10,
        alignment=TA_CENTER,
        fontName="Helvetica-Bold",
    )
)


# ============================================================
# Table helpers
# ============================================================

def paragraph_cell(
    value,
    style="TableText"
):

    return Paragraph(
        str(value),
        styles[style]
    )


def report_table(
    data,
    widths=None
):

    formatted = []

    for row_index, row in enumerate(data):

        formatted_row = []

        for value in row:

            style_name = (
                "TableHeader"
                if row_index == 0
                else "TableText"
            )

            formatted_row.append(
                paragraph_cell(
                    value,
                    style_name
                )
            )

        formatted.append(
            formatted_row
        )

    table = Table(
        formatted,
        colWidths=widths,
        repeatRows=1,
        hAlign="LEFT"
    )

    table.setStyle(
        TableStyle([
            (
                "BACKGROUND",
                (0, 0),
                (-1, 0),
                colors.HexColor(
                    "#E8EEF5"
                )
            ),
            (
                "GRID",
                (0, 0),
                (-1, -1),
                0.4,
                colors.HexColor(
                    "#BBBBBB"
                )
            ),
            (
                "VALIGN",
                (0, 0),
                (-1, -1),
                "TOP"
            ),
            (
                "LEFTPADDING",
                (0, 0),
                (-1, -1),
                6
            ),
            (
                "RIGHTPADDING",
                (0, 0),
                (-1, -1),
                6
            ),
            (
                "TOPPADDING",
                (0, 0),
                (-1, -1),
                5
            ),
            (
                "BOTTOMPADDING",
                (0, 0),
                (-1, -1),
                5
            ),
        ])
    )

    return table


def add_section(
    story,
    title
):

    story.append(
        Paragraph(
            title,
            styles["SectionTitle"]
        )
    )


def add_body(
    story,
    text
):

    story.append(
        Paragraph(
            text,
            styles["BodyReport"]
        )
    )


def add_note(
    story,
    text
):

    story.append(
        Paragraph(
            text,
            styles["Note"]
        )
    )


# ============================================================
# PDF document
# ============================================================

doc = SimpleDocTemplate(
    args.output_pdf,
    pagesize=A4,
    rightMargin=1.6 * cm,
    leftMargin=1.6 * cm,
    topMargin=1.5 * cm,
    bottomMargin=1.5 * cm,
)


story = []


# ============================================================
# Title
# ============================================================

story.append(
    Paragraph(
        "Metagenomic Analysis Report",
        styles["ReportTitle"]
    )
)

story.append(
    Paragraph(
        "Integrated Oxford Nanopore metagenomic analysis",
        styles["ReportSubtitle"]
    )
)

story.append(
    report_table(
        [
            [
                "Sample",
                "Result"
            ],
            [
                "Sample name",
                args.sample_id
            ]
        ],
        [
            5.0 * cm,
            10.0 * cm
        ]
    )
)

story.append(
    Spacer(
        1,
        0.3 * cm
    )
)

add_body(
    story,
    (
        "This report summarizes sequencing quality, assembly, "
        "gene prediction, genome reconstruction, taxonomic "
        "composition and functional potential identified from "
        "the analyzed metagenomic sample. Read-based and "
        "genome-resolved analyses are presented separately "
        "because they describe different aspects of the microbial "
        "community."
    )
)


# ============================================================
# Analysis overview
# ============================================================

add_section(
    story,
    "How to Understand This Analysis"
)

add_body(
    story,
    (
        "Oxford Nanopore sequencing generates individual DNA "
        "reads originating from organisms present in the sample. "
        "After quality filtering, overlapping reads are assembled "
        "into longer contiguous sequences. Protein-coding genes "
        "are predicted from these assembled sequences, while "
        "genome binning attempts to reconstruct genome-like "
        "groups of contigs."
    )
)

add_body(
    story,
    (
        "Taxonomy is examined from two complementary perspectives. "
        "Kraken2 and Bracken classify sequencing reads directly "
        "and therefore describe the overall read-based community. "
        "GTDB-Tk classifies reconstructed genome bins and CoverM "
        "estimates the abundance associated with those bins. "
        "Consequently, these two taxonomic views are not expected "
        "to produce identical abundance profiles."
    )
)


# ============================================================
# QC
# ============================================================

add_section(
    story,
    "Sequencing Quality Control"
)

qc_table = [
    [
        "Metric",
        "Raw reads",
        "Filtered reads"
    ],
    [
        "Number of reads",
        f"{raw_reads:,}",
        f"{filtered_reads:,}"
    ],
    [
        "Total bases",
        human_number(
            raw_bases
        ),
        human_number(
            filtered_bases
        )
    ],
    [
        "Mean read length",
        (
            f"{raw_qc.get('mean_length', 0):,.1f} bp"
        ),
        (
            f"{filtered_qc.get('mean_length', 0):,.1f} bp"
        )
    ],
    [
        "Read N50",
        (
            f"{raw_qc.get('read_length_n50', 0):,.0f} bp"
        ),
        (
            f"{filtered_qc.get('read_length_n50', 0):,.0f} bp"
        )
    ],
    [
        "Mean read quality",
        (
            f"{raw_qc.get('mean_quality', 0):.1f}"
        ),
        (
            f"{filtered_qc.get('mean_quality', 0):.1f}"
        )
    ],
]

story.append(
    report_table(
        qc_table,
        [
            5.0 * cm,
            5.0 * cm,
            5.0 * cm
        ]
    )
)

story.append(
    Spacer(
        1,
        0.25 * cm
    )
)

if os.path.exists(
    qc_figure
):

    story.append(
        Image(
            qc_figure,
            width=13.5 * cm,
            height=8.1 * cm
        )
    )

retained_reads_pct = (
    filtered_reads / raw_reads * 100
    if raw_reads
    else 0
)

retained_bases_pct = (
    filtered_bases / raw_bases * 100
    if raw_bases
    else 0
)

add_note(
    story,
    (
        f"After filtering, {retained_reads_pct:.1f}% of reads and "
        f"{retained_bases_pct:.1f}% of sequenced bases were retained. "
        "The difference between read retention and base retention "
        "indicates that filtering preferentially retained longer "
        "and/or higher-quality reads."
    )
)


# ============================================================
# Assembly
# ============================================================

add_section(
    story,
    "Assembly"
)

assembly_table = [
    [
        "Metric",
        "Result"
    ],
    [
        "Assembly sequences",
        assembly_stats.get(
            "num_seqs",
            ""
        )
    ],
    [
        "Total assembly length",
        assembly_stats.get(
            "sum_len",
            ""
        )
    ],
    [
        "Assembly N50",
        assembly_stats.get(
            "N50",
            ""
        )
    ],
    [
        "Largest sequence",
        assembly_stats.get(
            "max_len",
            ""
        )
    ],
]

story.append(
    report_table(
        assembly_table,
        [
            8.0 * cm,
            7.0 * cm
        ]
    )
)

add_note(
    story,
    (
        "Assembly combines overlapping sequencing reads into longer "
        "contiguous sequences. N50 summarizes assembly continuity: "
        "half of the assembled sequence length is contained in "
        "contigs at least as long as the reported N50 value."
    )
)


# ============================================================
# Gene catalog
# ============================================================

add_section(
    story,
    "Gene Prediction and Protein Catalog"
)

add_body(
    story,
    (
        "Protein-coding genes were predicted from the assembled "
        "metagenome. Similar proteins were subsequently clustered "
        "to reduce redundancy and produce a representative protein "
        "catalog for downstream functional analysis."
    )
)

gene_table = [
    [
        "Metric",
        "Result"
    ],
    [
        "Predicted proteins",
        protein_stats.get(
            "num_seqs",
            ""
        )
    ],
    [
        "Representative proteins",
        catalog_stats.get(
            "num_seqs",
            ""
        )
    ],
]

story.append(
    report_table(
        gene_table,
        [
            8.0 * cm,
            7.0 * cm
        ]
    )
)


# ============================================================
# Genome binning
# ============================================================

add_section(
    story,
    "Genome Binning"
)

add_body(
    story,
    (
        f"The assembly produced <b>{bin_count}</b> candidate genome "
        "bins. A bin represents a collection of assembled contigs "
        "predicted to originate from the same or a closely related "
        "microbial population. The number of bins should therefore "
        "not be interpreted directly as the number of species, "
        "because bins can be incomplete and multiple bins may "
        "represent related organisms."
    )
)


# ============================================================
# Read taxonomy
# ============================================================

story.append(
    PageBreak()
)

add_section(
    story,
    "Read-Based Community Composition"
)

add_body(
    story,
    (
        "Kraken2 classified the filtered sequencing reads against "
        "the taxonomic reference database, and Bracken refined "
        "these classifications to estimate species-level abundance. "
        "This provides the broadest taxonomic view of the microbial "
        "community because it does not require organisms to be "
        "successfully assembled into genome bins."
    )
)

if bracken_top:

    read_tax_table = [
        [
            "Species",
            "Estimated abundance (%)"
        ]
    ]

    for name, value in bracken_top[:10]:

        read_tax_table.append([
            name,
            f"{value:.2f}"
        ])

    story.append(
        report_table(
            read_tax_table,
            [
                11.0 * cm,
                4.0 * cm
            ]
        )
    )

    story.append(
        Spacer(
            1,
            0.25 * cm
        )
    )

    if os.path.exists(
        taxonomy_figure
    ):

        story.append(
            Image(
                taxonomy_figure,
                width=13.8 * cm,
                height=10.0 * cm
            )
        )

    if bracken_top:

        dominant_name = bracken_top[0][0]
        dominant_value = bracken_top[0][1]

        add_note(
            story,
            (
                f"The most abundant species-level assignment in "
                f"the read-based profile was <b>{dominant_name}</b> "
                f"at approximately <b>{dominant_value:.2f}%</b>. "
                "This value represents the proportion estimated "
                "from classified sequencing reads and should not "
                "be interpreted as an absolute cell count."
            )
        )


# ============================================================
# Genome-resolved taxonomy
# ============================================================

story.append(
    PageBreak()
)

add_section(
    story,
    "Genome-Resolved Taxonomy"
)

add_body(
    story,
    (
        "GTDB-Tk assigns taxonomy to reconstructed genome bins, "
        "while CoverM estimates the sequencing abundance associated "
        "with those bins. This analysis therefore describes only "
        "the portion of the community that was sufficiently "
        "assembled and grouped into genome-like populations."
    )
)

if gtdb_top:

    gtdb_table = [
        [
            "Reconstructed taxon",
            "Bins",
            "Relative abundance (%)"
        ]
    ]

    for species, abundance in gtdb_top[:10]:

        gtdb_table.append([
            species,
            species_bins.get(
                species,
                0
            ),
            f"{abundance:.2f}"
        ])

    story.append(
        report_table(
            gtdb_table,
            [
                9.0 * cm,
                2.0 * cm,
                4.0 * cm
            ]
        )
    )

    story.append(
        Spacer(
            1,
            0.25 * cm
        )
    )

    if os.path.exists(
        gtdb_figure
    ):

        story.append(
            Image(
                gtdb_figure,
                width=13.8 * cm,
                height=10.0 * cm
            )
        )

add_note(
    story,
    (
        "<b>Why can this differ from Kraken2/Bracken?</b> "
        "Read-based taxonomy examines individual sequencing reads, "
        "whereas genome-resolved taxonomy requires reads to first "
        "assemble into contigs and then form sufficiently coherent "
        "genome bins. Organisms that are abundant in the read-level "
        "analysis can therefore be absent, fragmented or represented "
        "differently among reconstructed bins."
    )
)


# ============================================================
# Functional prediction
# ============================================================

story.append(
    PageBreak()
)

add_section(
    story,
    "Functional Prediction"
)

add_body(
    story,
    (
        "KofamScan compares representative proteins with KEGG "
        "Orthology profile hidden Markov models. A threshold-passing "
        "assignment indicates that a protein showed sufficient "
        "evidence for assignment to a KEGG Orthology functional "
        "group."
    )
)

functional_overview = [
    [
        "Result",
        "Value"
    ],
    [
        "Representative proteins screened",
        f"{safe_int(catalog_stats.get('num_seqs', 0)):,}"
    ],
    [
        "Accepted KO assignments",
        f"{kofam_accepted_count:,}"
    ],
    [
        "Unique KO identifiers",
        f"{len(kofam_accepted_kos):,}"
    ],
]

story.append(
    report_table(
        functional_overview,
        [
            10.5 * cm,
            4.5 * cm
        ]
    )
)

story.append(
    Spacer(
        1,
        0.3 * cm
    )
)

if os.path.exists(
    kofam_figure
):

    story.append(
        Image(
            kofam_figure,
            width=14.5 * cm,
            height=9.8 * cm
        )
    )

add_note(
    story,
    (
        "The figure summarizes the most frequently detected "
        "threshold-passing functional assignments. The bar length "
        "represents how many accepted protein assignments were "
        "associated with each function. These results describe "
        "the predicted functional potential encoded by the "
        "metagenomic protein catalog; they do not demonstrate that "
        "the corresponding genes were actively expressed."
    )
)


# ============================================================
# Functional annotation
# ============================================================

story.append(
    PageBreak()
)

add_section(
    story,
    "Functional Annotation"
)

add_body(
    story,
    (
        "eggNOG-mapper provides complementary orthology-based "
        "functional annotation of the representative protein "
        "catalog. The results can be summarized at broader "
        "functional-category and metabolic-pathway levels, making "
        "the overall functional profile easier to interpret than "
        "individual annotation identifiers alone."
    )
)

annotation_overview = [
    [
        "Annotation coverage",
        "Proteins"
    ],
    [
        "Proteins with eggNOG annotation",
        f"{eggnog_count:,}"
    ],
    [
        "Proteins with COG category",
        f"{eggnog_cog_count:,}"
    ],
    [
        "Proteins with KEGG Orthology annotation",
        f"{eggnog_ko_count:,}"
    ],
    [
        "Proteins associated with KEGG pathway",
        f"{eggnog_pathway_count:,}"
    ],
]

story.append(
    report_table(
        annotation_overview,
        [
            10.5 * cm,
            4.5 * cm
        ]
    )
)

story.append(
    Spacer(
        1,
        0.3 * cm
    )
)


# ------------------------------------------------------------
# COG profile
# ------------------------------------------------------------

story.append(
    Paragraph(
        "Broad Functional Profile",
        styles["SubsectionTitle"]
    )
)

add_body(
    story,
    (
        "COG categories group annotated proteins according to "
        "broad biological roles. The distribution below therefore "
        "shows which general functional classes are most strongly "
        "represented in the predicted protein catalog."
    )
)

if os.path.exists(
    cog_figure
):

    story.append(
        Image(
            cog_figure,
            width=14.5 * cm,
            height=10.0 * cm
        )
    )

if eggnog_cog_plot_counts:

    top_cog_code, top_cog_count = max(
        eggnog_cog_plot_counts.items(),
        key=lambda item: item[1]
    )

    top_cog_name = COG_LABELS.get(
        top_cog_code,
        top_cog_code
    )

    add_note(
        story,
        (
            f"The largest COG category was "
            f"<b>{top_cog_name}</b>, containing "
            f"<b>{top_cog_count:,}</b> assignments. "
            "A large unknown-function category indicates that "
            "many predicted proteins have recognizable orthologous "
            "relationships but still lack a sufficiently specific "
            "biological function."
        )
    )


# ------------------------------------------------------------
# KEGG pathways
# ------------------------------------------------------------

story.append(
    Paragraph(
        "Metabolic and Cellular Pathway Profile",
        styles["SubsectionTitle"]
    )
)

add_body(
    story,
    (
        "KEGG pathway associations summarize how annotated proteins "
        "are distributed among known metabolic and cellular systems. "
        "Because a single protein can participate in more than one "
        "pathway, pathway counts are overlapping annotations rather "
        "than mutually exclusive groups."
    )
)

if os.path.exists(
    pathway_figure
):

    story.append(
        Image(
            pathway_figure,
            width=14.5 * cm,
            height=10.0 * cm
        )
    )

if eggnog_pathway_counts:

    top_pathway_code, top_pathway_count = (
        eggnog_pathway_counts.most_common(1)[0]
    )

    top_pathway_name = KEGG_PATHWAY_LABELS.get(
        top_pathway_code,
        top_pathway_code
    )

    add_note(
        story,
        (
            f"The most frequently represented pathway grouping was "
            f"<b>{top_pathway_name}</b> with "
            f"<b>{top_pathway_count:,}</b> protein associations. "
            "The pathway profile indicates broad metabolic and "
            "cellular capabilities encoded by the community. "
            "It should not be interpreted as direct evidence of "
            "pathway activity, which would require expression- or "
            "activity-based measurements."
        )
    )


# ============================================================
# Statistical / ecological analyses
# ============================================================

add_section(
    story,
    "Statistical and Ecological Analysis"
)

add_body(
    story,
    (
        "The current report represents one biological sample. "
        "Comparative analyses such as differential abundance, "
        "between-sample beta diversity, PERMANOVA and ecological "
        "association-network inference are therefore not performed. "
        "These analyses require an appropriate multi-sample design "
        "and, for inferential comparisons, biological replication."
    )
)


# ============================================================
# Overall interpretation
# ============================================================

add_section(
    story,
    "Overall Interpretation"
)

add_body(
    story,
    (
        "The sequencing dataset retained substantial sequence yield "
        "after quality filtering and produced a large metagenomic "
        "assembly containing predicted protein-coding genes. "
        "Genome binning reconstructed candidate microbial genome "
        "populations, while read-based taxonomy provided a broader "
        "view of community composition."
    )
)

add_body(
    story,
    (
        "The taxonomic results from Kraken2/Bracken and GTDB-Tk/"
        "CoverM should be interpreted as complementary rather than "
        "competing measurements. Kraken2 and Bracken describe "
        "organisms represented among sequencing reads, whereas "
        "GTDB-Tk and CoverM describe organisms represented among "
        "successfully reconstructed genome bins."
    )
)

add_body(
    story,
    (
        "Functional analysis identified thousands of KEGG Orthology "
        "groups and extensive orthology-based annotations across the "
        "representative protein catalog. The COG and KEGG pathway "
        "profiles provide an overview of the biological functions "
        "and metabolic capabilities encoded by the microbial "
        "community. These findings represent genomic functional "
        "potential and do not by themselves establish active gene "
        "expression or metabolic activity."
    )
)


# ============================================================
# Build PDF
# ============================================================

doc.build(
    story
)

print(
    f"Report written to {args.output_pdf}"
)

print(
    f"Summary written to {args.output_summary}"
)

