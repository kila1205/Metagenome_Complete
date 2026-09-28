# Complete Long-Read Metagenome Analysis Pipeline

A modular Nextflow pipeline for end-to-end analysis of Oxford Nanopore long-read metagenomic sequencing data.

The workflow integrates quality control, taxonomic profiling, metagenome assembly, gene prediction, genome binning, abundance estimation, protein clustering, genome-resolved taxonomy, functional prediction, functional annotation, and automatic report generation.

## Workflow Overview

ONT FASTQ  
→ Quality Control: NanoPlot → Filtlong → NanoPlot  
→ Read-Based Taxonomy: Kraken2 + Bracken  
→ Assembly: Flye + QUAST  
→ Gene Prediction: Prodigal  
→ Genome Binning: MetaBAT2  
→ Bin Abundance: CoverM  
→ Protein Clustering: MMseqs2  
→ Genome-Resolved Taxonomy: GTDB-Tk  
→ Functional Prediction: KofamScan  
→ Functional Annotation: eggNOG-mapper  
→ Automatic PDF Report

## Main Analyses

### Quality Control

Raw Oxford Nanopore reads are evaluated with NanoPlot, filtered using Filtlong, and evaluated again after filtering.

The QC stage provides read count, sequencing yield, read-length statistics, read N50, read quality, and GC content.

### Read-Based Taxonomy

Kraken2 classifies sequencing reads against a user-provided reference database. Bracken estimates taxonomic abundance from the Kraken2 classifications.

This analysis describes organisms represented directly among the sequencing reads.

### Metagenome Assembly

Flye performs long-read metagenome assembly. QUAST and SeqKit summarize assembly quality, including total assembly length, number of sequences, and N50.

### Gene Prediction

Prodigal predicts protein-coding genes from the assembled metagenomic sequences.

### Genome Binning

MetaBAT2 reconstructs candidate metagenome-assembled genome bins using sequence composition and coverage information.

### Abundance Estimation

CoverM estimates the relative abundance associated with reconstructed genome bins.

### Protein Clustering

MMseqs2 clusters predicted proteins to produce a non-redundant representative protein catalog for downstream functional analysis.

### Genome-Resolved Taxonomy

GTDB-Tk assigns taxonomy to reconstructed genome bins.

Kraken2/Bracken and GTDB-Tk answer different biological questions. Kraken2 and Bracken classify sequencing reads, whereas GTDB-Tk classifies genome bins that were successfully reconstructed from the assembly.

Therefore, their taxonomic profiles are not expected to be identical.

### Functional Prediction

KofamScan compares representative proteins against KEGG Orthology profile HMMs.

The analysis distinguishes candidate matches from threshold-passing accepted KO assignments.

### Functional Annotation

eggNOG-mapper provides large-scale orthology-based functional annotation, including COG categories, KEGG Orthology identifiers, and KEGG pathway associations.

Environmental metagenomes commonly contain proteins without well-characterized reference homologs. Unknown or incompletely annotated proteins therefore do not automatically indicate sequencing or pipeline failure.

## Automatic Report

The pipeline generates an integrated PDF report using:

`scripts/report/generate_metagenome_report.py`

The report summarizes:

- sequencing quality control
- assembly statistics
- predicted proteins
- genome binning
- read-based taxonomy
- genome-resolved taxonomy
- abundance
- KofamScan functional prediction
- eggNOG functional annotation
- COG functional profiles
- KEGG pathway profiles
- overall biological interpretation

A machine-readable TSV summary is generated together with the PDF report.

## Required Reference Databases

Large reference databases are not stored in this repository.

The following parameters must point to the appropriate local databases:

- `--kraken_db` — Kraken2 database
- `--gtdbtk_data_path` — GTDB-Tk database
- `--kofam_profiles` — KofamScan profile directory
- `--kofam_ko_list` — KofamScan KO list
- `--eggnog_data_dir` — eggNOG-mapper database

## Repository Structure

    Metagenome_Complete/
    ├── main.nf
    ├── nextflow.config
    ├── README.md
    ├── conf/
    ├── envs/
    ├── modules/
    ├── samplesheets/
    └── scripts/

Raw sequencing data, biological databases, Conda caches, Nextflow work directories, and analysis results are intentionally excluded from the repository.

## Example Run

For long HPC runs, using a persistent terminal session such as tmux is recommended.

Start tmux:

    tmux new -s metagenome

Example Nextflow execution:

    nextflow run main.nf \
        --input /path/to/sample.fastq.gz \
        --sample_id Sample01 \
        --outdir results \
        --kraken_db /path/to/kraken2/database \
        --gtdbtk_data_path /path/to/gtdbtk/database \
        --kofam_profiles /path/to/kofam/profiles \
        --kofam_ko_list /path/to/kofam/ko_list \
        --eggnog_data_dir /path/to/eggnog/database

## Major Software

The pipeline uses modular Conda environments under `envs/`.

Major tools include NanoPlot, Filtlong, MultiQC, Kraken2, Bracken, Flye, QUAST, SeqKit, Prodigal, minimap2, SAMtools, MetaBAT2, CoverM, MMseqs2, GTDB-Tk, KofamScan, eggNOG-mapper, ReportLab, and Matplotlib.

## Output

Major pipeline outputs include:

- QC reports and filtered reads
- metagenome assembly and assembly statistics
- predicted proteins
- candidate genome bins
- read-based taxonomic profiles
- reconstructed-bin abundance estimates
- GTDB-Tk genome taxonomy
- representative protein catalog
- KofamScan functional predictions
- eggNOG functional annotations
- automatic PDF report
- summary TSV
- Nextflow execution information

## Interpretation Notes

Read-based taxonomy and genome-resolved taxonomy are complementary analyses and should not be expected to produce identical abundance profiles.

A high proportion of unclassified reads, unknown proteins, or incompletely classified genome bins does not automatically indicate poor sequencing quality. Environmental microbial communities often contain organisms and genes that are incompletely represented in current reference databases.

Functional annotations represent the portion of the predicted protein catalog that can currently be associated with known reference functions.

## HPC Configuration

`conf/cluster.config` contains settings intended for HPC execution.

Cluster-specific paths, resource limits, scheduler configuration, Conda configuration, and database locations may need to be adapted when running the pipeline on another system.

## Status

This pipeline has been developed and tested for long-read environmental metagenomic analysis using Oxford Nanopore sequencing data.
