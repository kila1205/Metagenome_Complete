#!/usr/bin/env python3

from pathlib import Path
import re
import csv
from datetime import datetime

PROJECT = Path("/media/bio-pool/KIL/Metagenome/Metagenome_Complete_GitHub")
OUTDIR = PROJECT / "resource_report"
OUTDIR.mkdir(parents=True, exist_ok=True)

TXT = OUTDIR / "metagenome_compute_resources.txt"
TSV = OUTDIR / "metagenome_compute_resources.tsv"

# Friendly descriptions for known processes.
DESCRIPTIONS = {
    "RAW_NANOPLOT":       ("Raw-read QC", "NanoPlot"),
    "FILTLONG":           ("Read filtering", "Filtlong"),
    "FILTERED_NANOPLOT":  ("Filtered-read QC", "NanoPlot"),
    "MULTIQC":            ("QC summary", "MultiQC"),
    "MERGE_READS":        ("Read merging", "pigz / shell"),
    "METAFLYE":           ("Metagenome assembly", "metaFlye"),
    "ASSEMBLY_STATS":     ("Assembly statistics", "Assembly statistics"),
    "QUAST":              ("Assembly quality assessment", "QUAST"),
    "PRODIGAL":           ("Gene prediction", "Prodigal"),
    "PROTEIN_STATS":      ("Protein statistics", "Protein statistics"),
    "MAP_READS":          ("Read mapping for binning", "Minimap2 / Samtools"),
    "CONTIG_DEPTH":       ("Contig depth calculation", "MetaBAT2 utilities"),
    "METABAT2":           ("Genome binning", "MetaBAT2"),
    "BIN_STATS":          ("Bin statistics", "Bin statistics"),
    "KRAKEN2_CLASSIFY":   ("Read-based taxonomy", "Kraken2"),
    "BRACKEN_ABUNDANCE":  ("Taxonomic abundance estimation", "Bracken"),
    "COVERM_ABUNDANCE":   ("Genome/bin abundance", "CoverM"),
    "COMBINE_PROTEINS":   ("Protein catalog preparation", "Shell"),
    "MMSEQS_CLUSTER":     ("Protein clustering", "MMseqs2"),
    "CATALOG_STATS":      ("Protein catalog statistics", "Catalog statistics"),
    "GTDBTK_CLASSIFY":    ("Genome taxonomy", "GTDB-Tk"),
    "PARSE_TAXONOMY":     ("Taxonomy result parsing", "Python"),
    "KOFAMSCAN":          ("Functional prediction", "KofamScan"),
    "EGGNOG_MAPPER":      ("Functional annotation", "eggNOG-mapper"),
    "GENERATE_REPORT":    ("Final report generation", "Report generator"),
}

# Conceptual biological order.
ORDER = [
    "RAW_NANOPLOT",
    "FILTLONG",
    "FILTERED_NANOPLOT",
    "MULTIQC",
    "MERGE_READS",
    "METAFLYE",
    "ASSEMBLY_STATS",
    "QUAST",
    "PRODIGAL",
    "PROTEIN_STATS",
    "MAP_READS",
    "CONTIG_DEPTH",
    "METABAT2",
    "BIN_STATS",
    "KRAKEN2_CLASSIFY",
    "BRACKEN_ABUNDANCE",
    "COVERM_ABUNDANCE",
    "COMBINE_PROTEINS",
    "MMSEQS_CLUSTER",
    "CATALOG_STATS",
    "GTDBTK_CLASSIFY",
    "PARSE_TAXONOMY",
    "KOFAMSCAN",
    "EGGNOG_MAPPER",
    "GENERATE_REPORT",
]

def extract_block(text, start):
    brace = text.find("{", start)
    if brace == -1:
        return None

    depth = 0

    for i in range(brace, len(text)):
        if text[i] == "{":
            depth += 1
        elif text[i] == "}":
            depth -= 1
            if depth == 0:
                return text[start:i+1]

    return None


def clean(v):
    if v is None:
        return None

    v = v.split("//")[0].strip().rstrip(";").strip()

    if len(v) >= 2 and (
        (v[0] == "'" and v[-1] == "'") or
        (v[0] == '"' and v[-1] == '"')
    ):
        v = v[1:-1]

    return v


def directive(block, name):
    patterns = [
        rf"(?m)^\s*{name}\s*=\s*(.+?)\s*$",
        rf"(?m)^\s*{name}\s+(.+?)\s*$",
    ]

    for pattern in patterns:
        m = re.search(pattern, block)
        if m:
            return clean(m.group(1))

    return None


# ------------------------------------------------------------
# Discover all process declarations in modules
# ------------------------------------------------------------

processes = {}

for path in sorted((PROJECT / "modules").rglob("*.nf")):

    text = path.read_text(errors="replace")

    for match in re.finditer(
        r"\bprocess\s+([A-Za-z0-9_]+)\s*\{",
        text
    ):
        name = match.group(1)

        block = extract_block(text, match.start())

        if not block:
            continue

        processes[name] = {
            "process": name,
            "source": str(path.relative_to(PROJECT)),
            "cpus": directive(block, "cpus"),
            "memory": directive(block, "memory"),
            "time": directive(block, "time"),
        }


# ------------------------------------------------------------
# Read withName overrides from configuration files.
# These values override module-level values.
# ------------------------------------------------------------

config_paths = [
    PROJECT / "nextflow.config",
    PROJECT / "conf" / "cluster.config",
]

shared_config = Path(
    "/media/bio-pool/software/nextflow/slurm-bersama.config"
)

if shared_config.exists():
    config_paths.append(shared_config)

for path in config_paths:

    if not path.exists():
        continue

    text = path.read_text(errors="replace")

    pattern = re.compile(
        r"withName\s*:\s*['\"]([^'\"]+)['\"]\s*\{"
    )

    for match in pattern.finditer(text):

        name = match.group(1)
        block = extract_block(text, match.start())

        if not block:
            continue

        if name not in processes:
            processes[name] = {
                "process": name,
                "source": "",
                "cpus": None,
                "memory": None,
                "time": None,
            }

        c = directive(block, "cpus")
        m = directive(block, "memory")
        t = directive(block, "time")

        if c is not None:
            processes[name]["cpus"] = c

        if m is not None:
            processes[name]["memory"] = m

        if t is not None:
            processes[name]["time"] = t

        rel = (
            str(path.relative_to(PROJECT))
            if PROJECT in path.parents
            else str(path)
        )

        if c or m or t:
            processes[name]["source"] += (
                ("; " if processes[name]["source"] else "")
                + rel
            )


# ------------------------------------------------------------
# Global defaults
# ------------------------------------------------------------

global_defaults = {
    "cpus": None,
    "memory": None,
    "time": None,
}

nextflow_config = PROJECT / "nextflow.config"

if nextflow_config.exists():

    text = nextflow_config.read_text(errors="replace")

    # Locate broad process block if present.
    for match in re.finditer(r"\bprocess\s*\{", text):

        block = extract_block(text, match.start())

        if not block:
            continue

        for key in global_defaults:
            value = directive(block, key)

            if value is not None:
                global_defaults[key] = value


# ------------------------------------------------------------
# Global resource limits
# ------------------------------------------------------------

resource_limits = "Not explicitly defined"

cluster = PROJECT / "conf" / "cluster.config"

if cluster.exists():

    text = cluster.read_text(errors="replace")

    m = re.search(
        r"resourceLimits\s*=\s*\[(.*?)\]",
        text,
        flags=re.S
    )

    if m:
        resource_limits = " ".join(m.group(1).split())


# ------------------------------------------------------------
# Prepare report rows
# ------------------------------------------------------------

rows = []

ordered_names = []

for name in ORDER:
    if name in processes:
        ordered_names.append(name)

for name in sorted(processes):
    if name not in ordered_names:
        ordered_names.append(name)


for name in ordered_names:

    p = processes[name]

    analysis, tool = DESCRIPTIONS.get(
        name,
        ("Pipeline process", name)
    )

    cpus = p["cpus"]
    memory = p["memory"]
    walltime = p["time"]

    cpu_source = "Process-specific"
    mem_source = "Process-specific"
    time_source = "Process-specific"

    if cpus is None:
        cpus = global_defaults["cpus"]
        cpu_source = "Global/default"

    if memory is None:
        memory = global_defaults["memory"]
        mem_source = "Global/default"

    if walltime is None:
        walltime = global_defaults["time"]
        time_source = "Global/default"

    if cpus is None:
        cpus = "Nextflow/cluster default"
        cpu_source = "No explicit value"

    if memory is None:
        memory = "Nextflow/cluster default"
        mem_source = "No explicit value"

    if walltime is None:
        walltime = "Scheduler/default"
        time_source = "No explicit value"

    rows.append({
        "process": name,
        "analysis": analysis,
        "tool": tool,
        "cpus": cpus,
        "memory": memory,
        "time": walltime,
        "cpu_source": cpu_source,
        "memory_source": mem_source,
        "time_source": time_source,
        "source": p["source"] or "Module/default configuration",
    })


# ------------------------------------------------------------
# TSV
# ------------------------------------------------------------

with TSV.open("w", newline="") as fh:

    writer = csv.writer(fh, delimiter="\t")

    writer.writerow([
        "Pipeline Process",
        "Analysis",
        "Tool",
        "CPU",
        "Memory",
        "Walltime",
        "CPU Source",
        "Memory Source",
        "Walltime Source",
        "Configuration Source",
    ])

    for r in rows:
        writer.writerow([
            r["process"],
            r["analysis"],
            r["tool"],
            r["cpus"],
            r["memory"],
            r["time"],
            r["cpu_source"],
            r["memory_source"],
            r["time_source"],
            r["source"],
        ])


# ------------------------------------------------------------
# TXT report
# ------------------------------------------------------------

with TXT.open("w") as out:

    out.write("=" * 88 + "\n")
    out.write("METAGENOME BIOINFORMATICS PIPELINE\n")
    out.write("COMPUTATIONAL RESOURCE REQUIREMENT REPORT\n")
    out.write("=" * 88 + "\n\n")

    out.write(
        f"Generated : "
        f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n"
    )

    out.write(f"Pipeline  : {PROJECT}\n")
    out.write("Workflow  : Nextflow DSL2\n")
    out.write("Executor  : Slurm\n")
    out.write("Queue     : long\n")
    out.write("QoS       : normal\n\n")

    out.write("PURPOSE\n")
    out.write("-" * 88 + "\n")

    out.write(
        "This report documents the computational resource configuration of "
        "the metagenome\n"
        "bioinformatics pipeline. CPU, memory and walltime values are "
        "extracted directly\n"
        "from the Nextflow modules and configuration files. These values "
        "represent configured\n"
        "resource requests rather than measured peak hardware consumption.\n\n"
    )

    out.write("PROCESS RESOURCE REQUIREMENTS\n")
    out.write("=" * 88 + "\n\n")

    for r in rows:

        out.write(f"Process       : {r['process']}\n")
        out.write(f"Analysis      : {r['analysis']}\n")
        out.write(f"Tool          : {r['tool']}\n")
        out.write(f"CPU           : {r['cpus']}\n")
        out.write(f"Memory        : {r['memory']}\n")
        out.write(f"Walltime      : {r['time']}\n")
        out.write(f"Source        : {r['source']}\n")
        out.write("-" * 88 + "\n")


    out.write("\nRESOURCE INTERPRETATION\n")
    out.write("=" * 88 + "\n")

    out.write(
        "A value reported as 'Nextflow/cluster default' or "
        "'Scheduler/default' means that\n"
        "the process does not declare a dedicated value in the pipeline "
        "configuration. In\n"
        "that case, execution uses the applicable Nextflow, Slurm, or shared "
        "cluster default.\n"
        "This does NOT mean that the process requires zero CPU or zero memory.\n\n"
    )

    out.write("GLOBAL EXECUTION CONFIGURATION\n")
    out.write("=" * 88 + "\n")
    out.write("Workflow manager      : Nextflow\n")
    out.write("Executor              : Slurm\n")
    out.write("Queue                 : long\n")
    out.write("QoS                   : normal\n")
    out.write("Configured queue size : 2 concurrent submitted jobs per workflow\n")
    out.write(f"Resource limits       : {resource_limits}\n\n")

    out.write(
        "Resource allocation is process-specific. The maximum configured "
        "resource limit is\n"
        "therefore not continuously consumed by the pipeline. Individual "
        "processes request\n"
        "resources independently through Slurm.\n"
    )


print("=" * 72)
print("REPORT GENERATED")
print("=" * 72)
print(f"TXT : {TXT}")
print(f"TSV : {TSV}")
