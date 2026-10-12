#!/usr/bin/env python3

import argparse
import os
import csv

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
)


def parse_seqkit(path):
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        return next(reader)


def fmt_int(value):
    return f"{int(float(value)):,}"


def fmt_float(value, digits=1):
    return f"{float(value):,.{digits}f}"


def build_report(args):
    assembly = parse_seqkit(args.assembly_stats)
    protein = parse_seqkit(args.protein_stats)

    os.makedirs(os.path.dirname(args.output_pdf), exist_ok=True)

    styles = getSampleStyleSheet()

    title_style = ParagraphStyle(
        "TitleCustom",
        parent=styles["Title"],
        fontName="Helvetica-Bold",
        fontSize=18,
        leading=22,
        alignment=TA_CENTER,
        spaceAfter=12,
    )

    subtitle_style = ParagraphStyle(
        "Subtitle",
        parent=styles["Normal"],
        fontName="Helvetica",
        fontSize=10,
        leading=14,
        alignment=TA_CENTER,
        textColor=colors.HexColor("#555555"),
        spaceAfter=20,
    )

    heading = ParagraphStyle(
        "HeadingCustom",
        parent=styles["Heading2"],
        fontName="Helvetica-Bold",
        fontSize=13,
        leading=16,
        spaceBefore=12,
        spaceAfter=8,
    )

    body = ParagraphStyle(
        "BodyCustom",
        parent=styles["BodyText"],
        fontName="Helvetica",
        fontSize=10,
        leading=15,
        spaceAfter=8,
    )

    note = ParagraphStyle(
        "Note",
        parent=body,
        fontSize=9,
        leading=13,
        textColor=colors.HexColor("#555555"),
    )

    doc = SimpleDocTemplate(
        args.output_pdf,
        pagesize=A4,
        rightMargin=1.8 * cm,
        leftMargin=1.8 * cm,
        topMargin=1.8 * cm,
        bottomMargin=1.8 * cm,
    )

    story = []

    story.append(Paragraph(
        "Metagenome Analysis Report",
        title_style
    ))
    story.append(Paragraph(
        f"Sample: {args.sample_id}",
        subtitle_style
    ))

    story.append(Paragraph("Analysis Overview", heading))

    story.append(Paragraph(
        "This report summarizes the metagenomic analysis results that were "
        "successfully generated for this sample. Quality-control, assembly, "
        "gene-prediction and other independent workflow stages were executed. "
        "Genome binning was also attempted; however, no genome bins were "
        "recovered from the assembled contigs. Consequently, analyses that "
        "require genome-bin FASTA files could not be completed.",
        body,
    ))

    overview_data = [
        ["Analysis component", "Status"],
        ["Sequence processing and quality control", "Completed"],
        ["Assembly", "Completed"],
        ["Gene prediction", "Completed"],
        ["Genome binning", "Completed - 0 bins recovered"],
        ["Bin abundance (CoverM)", "Not available - no genome bins"],
        ["Bin taxonomy (GTDB-Tk)", "Not available - no genome bins"],
    ]

    table = Table(overview_data, colWidths=[9.5 * cm, 6.2 * cm])
    table.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#E8E8E8")),
        ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
        ("FONTNAME", (0, 1), (-1, -1), "Helvetica"),
        ("FONTSIZE", (0, 0), (-1, -1), 9),
        ("LEADING", (0, 0), (-1, -1), 12),
        ("GRID", (0, 0), (-1, -1), 0.4, colors.HexColor("#AAAAAA")),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("LEFTPADDING", (0, 0), (-1, -1), 6),
        ("RIGHTPADDING", (0, 0), (-1, -1), 6),
        ("TOPPADDING", (0, 0), (-1, -1), 5),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 5),
    ]))

    story.append(table)
    story.append(Spacer(1, 12))

    story.append(Paragraph("Assembly Results", heading))

    assembly_data = [
        ["Metric", "Result"],
        ["Number of contigs", fmt_int(assembly["num_seqs"])],
        ["Total assembly length", f'{fmt_int(assembly["sum_len"])} bp'],
        ["Average contig length", f'{fmt_float(assembly["avg_len"])} bp'],
        ["Largest contig", f'{fmt_int(assembly["max_len"])} bp'],
        ["N50", f'{fmt_int(assembly["N50"])} bp'],
        ["GC content", f'{float(assembly["GC(%)"]):.2f}%'],
    ]

    table = Table(assembly_data, colWidths=[8 * cm, 7.7 * cm])
    table.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#E8E8E8")),
        ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
        ("FONTNAME", (0, 1), (-1, -1), "Helvetica"),
        ("FONTSIZE", (0, 0), (-1, -1), 9),
        ("GRID", (0, 0), (-1, -1), 0.4, colors.HexColor("#AAAAAA")),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("LEFTPADDING", (0, 0), (-1, -1), 6),
        ("TOPPADDING", (0, 0), (-1, -1), 5),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 5),
    ]))

    story.append(table)

    story.append(Paragraph(
        "The assembly process completed successfully and generated assembled "
        "contigs that could be used for subsequent gene prediction and "
        "binning. Assembly completion does not by itself indicate recovery "
        "of a complete microbial genome.",
        note,
    ))

    story.append(Paragraph("Gene Prediction", heading))

    protein_data = [
        ["Metric", "Result"],
        ["Predicted protein sequences", fmt_int(protein["num_seqs"])],
        ["Total protein residues", fmt_int(protein["sum_len"])],
        ["Average protein length", f'{fmt_float(protein["avg_len"])} aa'],
        ["Longest predicted protein", f'{fmt_int(protein["max_len"])} aa'],
    ]

    table = Table(protein_data, colWidths=[8 * cm, 7.7 * cm])
    table.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#E8E8E8")),
        ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
        ("FONTNAME", (0, 1), (-1, -1), "Helvetica"),
        ("FONTSIZE", (0, 0), (-1, -1), 9),
        ("GRID", (0, 0), (-1, -1), 0.4, colors.HexColor("#AAAAAA")),
        ("LEFTPADDING", (0, 0), (-1, -1), 6),
        ("TOPPADDING", (0, 0), (-1, -1), 5),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 5),
    ]))

    story.append(table)

    story.append(Paragraph("Genome Binning", heading))

    story.append(Paragraph(
        "MetaBAT2 completed the genome-binning stage, but no genome bins "
        "were recovered from this sample. The bin statistics file therefore "
        "contains no genome-bin entries.",
        body,
    ))

    bin_data = [
        ["Metric", "Result"],
        ["Genome bins recovered", "0"],
        ["Bin abundance analysis", "Not available"],
        ["GTDB-Tk bin classification", "Not available"],
    ]

    table = Table(bin_data, colWidths=[8 * cm, 7.7 * cm])
    table.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#E8E8E8")),
        ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
        ("FONTNAME", (0, 1), (-1, -1), "Helvetica"),
        ("FONTSIZE", (0, 0), (-1, -1), 9),
        ("GRID", (0, 0), (-1, -1), 0.4, colors.HexColor("#AAAAAA")),
        ("LEFTPADDING", (0, 0), (-1, -1), 6),
        ("TOPPADDING", (0, 0), (-1, -1), 5),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 5),
    ]))

    story.append(table)

    story.append(Paragraph("Interpretation", heading))

    story.append(Paragraph(
        "The sample generated an assembly and predicted protein sequences, "
        "showing that the sequence data were processable by the workflow. "
        "However, the available assembly did not support recovery of genome "
        "bins by MetaBAT2. Because CoverM bin-abundance analysis and GTDB-Tk "
        "classification depend on recovered genome bins, these results are "
        "not reported for this sample.",
        body,
    ))

    story.append(Paragraph(
        "The absence of recovered bins should not be interpreted as evidence "
        "that microorganisms were absent from the sample. It indicates that "
        "the available sequence and assembly information was insufficient "
        "for this workflow to reconstruct usable genome bins.",
        body,
    ))

    story.append(Paragraph("Report Scope", heading))

    story.append(Paragraph(
        "This is a partial analysis report. Only results produced by "
        "successfully completed and interpretable workflow stages are "
        "included. Missing bin-dependent results were not estimated or "
        "reconstructed.",
        body,
    ))

    doc.build(story)

    with open(args.output_summary, "w") as out:
        out.write("sample_id\tassembly_contigs\tassembly_length_bp\tN50_bp\t"
                  "largest_contig_bp\tGC_percent\tpredicted_proteins\tgenome_bins\tstatus\n")

        out.write(
            f'{args.sample_id}\t'
            f'{assembly["num_seqs"]}\t'
            f'{assembly["sum_len"]}\t'
            f'{assembly["N50"]}\t'
            f'{assembly["max_len"]}\t'
            f'{assembly["GC(%)"]}\t'
            f'{protein["num_seqs"]}\t'
            f'0\t'
            f'PARTIAL_NO_GENOME_BINS\n'
        )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--assembly-stats", required=True)
    parser.add_argument("--protein-stats", required=True)
    parser.add_argument("--output-pdf", required=True)
    parser.add_argument("--output-summary", required=True)

    args = parser.parse_args()
    build_report(args)


if __name__ == "__main__":
    main()
