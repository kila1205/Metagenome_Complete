#!/usr/bin/env python3

import argparse
import csv
import os
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER, TA_LEFT, TA_JUSTIFY
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import cm
from reportlab.platypus import (
    SimpleDocTemplate,
    Paragraph,
    Spacer,
    Table,
    TableStyle,
    Image,
)


# ============================================================
# Arguments
# ============================================================

parser = argparse.ArgumentParser(
    description="Generate partial metagenomic analysis reports."
)

parser.add_argument("--sample-id", required=True)
parser.add_argument("--raw-nanostats", required=True)
parser.add_argument("--filtered-nanostats", required=True)
parser.add_argument("--assembly-stats")
parser.add_argument(
    "--stage",
    required=True,
    choices=["filtering", "assembly"]
)
parser.add_argument("--output-pdf", required=True)
parser.add_argument("--output-summary", required=True)
parser.add_argument("--figures-dir", required=True)

args = parser.parse_args()

os.makedirs(args.figures_dir, exist_ok=True)
os.makedirs(os.path.dirname(args.output_pdf), exist_ok=True)
os.makedirs(os.path.dirname(args.output_summary), exist_ok=True)


# ============================================================
# Helpers
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


def human_number(number):
    number = float(number)

    if number >= 1_000_000_000:
        return f"{number / 1_000_000_000:.2f} Gb"
    if number >= 1_000_000:
        return f"{number / 1_000_000:.2f} Mb"
    if number >= 1_000:
        return f"{number / 1_000:.2f} kb"

    return f"{number:.0f}"


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


def read_seqkit_stats(path):
    if not path or not os.path.exists(path):
        return {}

    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        lines = [line for line in handle if line.strip()]

    if not lines:
        return {}

    reader = csv.DictReader(lines, delimiter="\t")
    rows = list(reader)

    return rows[0] if rows else {}


# ============================================================
# Read results
# ============================================================

raw_qc = parse_nanostats(args.raw_nanostats)
filtered_qc = parse_nanostats(args.filtered_nanostats)

raw_reads = raw_qc.get("number_reads", 0)
raw_bases = raw_qc.get("total_bases", 0)

filtered_reads = filtered_qc.get("number_reads", 0)
filtered_bases = filtered_qc.get("total_bases", 0)

assembly_stats = {}

if args.stage == "assembly":
    assembly_stats = read_seqkit_stats(args.assembly_stats)


# ============================================================
# QC figure
# ============================================================

qc_figure = os.path.join(
    args.figures_dir,
    "qc_comparison.png"
)

fig, ax = plt.subplots(figsize=(7.2, 4.3))

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

x = range(len(labels))

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

ax.set_xticks(list(x))
ax.set_xticklabels(labels)
ax.set_ylabel("Amount of sequencing data")
ax.set_title("Sequencing data before and after filtering")
ax.legend()

fig.tight_layout()

fig.savefig(
    qc_figure,
    dpi=180,
    bbox_inches="tight"
)

plt.close(fig)


# ============================================================
# Summary TSV
# ============================================================

summary_rows = [
    ("sample_id", args.sample_id),
    ("report_status", "Partial analysis"),
    ("last_reported_stage", args.stage),
    ("raw_reads", raw_reads),
    ("raw_bases", raw_bases),
    ("filtered_reads", filtered_reads),
    ("filtered_bases", filtered_bases),
]

if args.stage == "assembly":
    summary_rows.extend([
        (
            "assembly_sequences",
            assembly_stats.get("num_seqs", "")
        ),
        (
            "assembly_length",
            assembly_stats.get("sum_len", "")
        ),
        (
            "assembly_n50",
            assembly_stats.get("N50", "")
        ),
        (
            "assembly_max_length",
            assembly_stats.get("max_len", "")
        ),
    ])

with open(
    args.output_summary,
    "w",
    encoding="utf-8",
    newline=""
) as handle:

    writer = csv.writer(handle, delimiter="\t")
    writer.writerow(["metric", "value"])
    writer.writerows(summary_rows)


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
# PDF helpers
# ============================================================

def paragraph_cell(value, style="TableText"):
    return Paragraph(
        str(value),
        styles[style]
    )


def report_table(data, widths=None):
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
                paragraph_cell(value, style_name)
            )

        formatted.append(formatted_row)

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
                colors.HexColor("#E8EEF5")
            ),
            (
                "GRID",
                (0, 0),
                (-1, -1),
                0.4,
                colors.HexColor("#BBBBBB")
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


def add_section(story, title):
    story.append(
        Paragraph(
            title,
            styles["SectionTitle"]
        )
    )


def add_body(story, text):
    story.append(
        Paragraph(
            text,
            styles["BodyReport"]
        )
    )


def add_note(story, text):
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
            ["Sample", "Result"],
            ["Sample name", args.sample_id],
            ["Analysis status", "Partial analysis"],
        ],
        [
            5.0 * cm,
            10.0 * cm
        ]
    )
)

story.append(Spacer(1, 0.3 * cm))

if args.stage == "filtering":
    add_body(
        story,
        (
            "This report summarizes the sequencing quality-control "
            "and read-filtering results obtained for the analyzed "
            "Oxford Nanopore sample. Assembly was not successfully "
            "obtained; therefore, assembly-dependent downstream "
            "analyses are not included in this report."
        )
    )

else:
    add_body(
        story,
        (
            "This report summarizes the sequencing quality-control, "
            "read-filtering and assembly results obtained for the "
            "analyzed Oxford Nanopore sample. The report is limited "
            "to the validated analysis stages requested for this "
            "sample."
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
        "Oxford Nanopore sequencing generates individual DNA reads "
        "from the biological material present in the sample. Quality "
        "control is used to characterize these reads, after which "
        "filtering retains sequencing data suitable for subsequent "
        "analysis."
    )
)

if args.stage == "assembly":
    add_body(
        story,
        (
            "During assembly, overlapping sequencing reads are "
            "combined into longer contiguous sequences called contigs. "
            "Assembly statistics describe the amount and continuity "
            "of reconstructed sequence obtained from the available "
            "sequencing data."
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
        human_number(raw_bases),
        human_number(filtered_bases)
    ],
    [
        "Mean read length",
        f"{raw_qc.get('mean_length', 0):,.1f} bp",
        f"{filtered_qc.get('mean_length', 0):,.1f} bp"
    ],
    [
        "Read N50",
        f"{raw_qc.get('read_length_n50', 0):,.0f} bp",
        f"{filtered_qc.get('read_length_n50', 0):,.0f} bp"
    ],
    [
        "Mean read quality",
        f"{raw_qc.get('mean_quality', 0):.1f}",
        f"{filtered_qc.get('mean_quality', 0):.1f}"
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

story.append(Spacer(1, 0.25 * cm))

if os.path.exists(qc_figure):
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
        "Filtering is intended to retain reads that are more suitable "
        "for downstream sequence analysis."
    )
)


# ============================================================
# Barcode07: stop after filtering
# ============================================================

if args.stage == "filtering":

    add_section(
        story,
        "Analysis Status"
    )

    add_body(
        story,
        (
            "Quality-control and read-filtering steps were completed "
            "successfully. After filtering, the dataset contained "
            f"<b>{filtered_reads:,}</b> reads representing "
            f"<b>{human_number(filtered_bases)}</b> of sequence data."
        )
    )

    add_body(
        story,
        (
            "Metagenomic assembly was attempted but did not complete "
            "successfully for this dataset. Because a valid assembly "
            "was not obtained, assembly-dependent analyses are not "
            "reported here."
        )
    )

    add_note(
        story,
        (
            "<b>Reporting scope:</b> This report intentionally ends "
            "at the quality-control and read-filtering stage. "
            "Unavailable downstream results have not been estimated "
            "or inferred."
        )
    )


# ============================================================
# Barcode08: include assembly
# ============================================================

if args.stage == "assembly":

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
            "Assembly combines overlapping sequencing reads into "
            "longer contiguous sequences. N50 summarizes assembly "
            "continuity: half of the assembled sequence length is "
            "contained in contigs at least as long as the reported "
            "N50 value."
        )
    )

    add_section(
        story,
        "Analysis Status"
    )

    add_body(
        story,
        (
            "Quality control, read filtering and metagenomic assembly "
            "were completed for this sample. The assembly contained "
            f"<b>{assembly_stats.get('num_seqs', '')}</b> sequences "
            "with a total assembled length of "
            f"<b>{assembly_stats.get('sum_len', '')} bp</b> and an "
            f"N50 of <b>{assembly_stats.get('N50', '')} bp</b>."
        )
    )

    add_body(
        story,
        (
            "Genome binning did not recover usable genome bins from "
            "this assembly. Consequently, bin-dependent results such "
            "as genome-resolved taxonomic classification and bin-level "
            "abundance estimation are not included in this report."
        )
    )

    add_note(
        story,
        (
            "<b>Reporting scope:</b> This report intentionally "
            "presents results through the assembly stage. Results "
            "requiring successfully reconstructed genome bins are "
            "not reported."
        )
    )


# ============================================================
# Build PDF
# ============================================================

doc.build(story)

print(f"PDF report written to: {args.output_pdf}")
print(f"Summary written to: {args.output_summary}")
