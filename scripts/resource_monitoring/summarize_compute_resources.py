#!/usr/bin/env python3

import csv
import sys
from pathlib import Path
from statistics import mean


if len(sys.argv) != 5:
    print(
        "Usage: summarize_compute_resources.py "
        "SAMPLE TRACE HARDWARE OUTPUT_DIR"
    )
    sys.exit(1)


sample = sys.argv[1]
trace_file = Path(sys.argv[2])
hardware_file = Path(sys.argv[3])
output_dir = Path(sys.argv[4])

output_dir.mkdir(parents=True, exist_ok=True)

process_output = output_dir / "per_process_resources.tsv"
summary_output = output_dir / "compute_resource_summary.tsv"


def numeric(values):
    result = []

    for value in values:
        if value in (None, "", "NA", "N/A"):
            continue

        try:
            result.append(float(value))
        except ValueError:
            pass

    return result


# ============================================================
# HARDWARE TELEMETRY
# ============================================================

hardware = []

if hardware_file.exists():
    with hardware_file.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        hardware.extend(reader)


cpu_temps = []

for field in ("CPU1_temp_C", "CPU2_temp_C"):
    cpu_temps.extend(
        numeric([row.get(field) for row in hardware])
    )


ddr_fields = [
    "DDR4_A_temp_C",
    "DDR4_B_temp_C",
    "DDR4_C_temp_C",
    "DDR4_D_temp_C",
    "DDR4_I_temp_C",
    "DDR4_J_temp_C",
    "DDR4_K_temp_C",
    "DDR4_L_temp_C",
]

ddr_temps = []

for field in ddr_fields:
    ddr_temps.extend(
        numeric([row.get(field) for row in hardware])
    )


ram_values = numeric(
    [row.get("ram_used_GB") for row in hardware]
)


gpu_temps = []

for field in (
    "GPU0_temp_C",
    "GPU1_temp_C",
    "GPU2_temp_C",
):
    gpu_temps.extend(
        numeric([row.get(field) for row in hardware])
    )


# ============================================================
# NEXTFLOW TRACE
# ============================================================

trace_rows = []

if trace_file.exists():
    with trace_file.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        trace_rows.extend(reader)


# ============================================================
# PER-PROCESS TABLE
# ============================================================

process_fields = [
    "sample",
    "process",
    "task_id",
    "native_id",
    "status",
    "exit",
    "submit",
    "start",
    "complete",
    "duration",
    "realtime",
    "cpus",
    "memory",
    "percent_cpu",
    "peak_rss",
    "peak_vmem",
    "read_bytes",
    "write_bytes",
]


with process_output.open("w", newline="") as fh:

    writer = csv.DictWriter(
        fh,
        fieldnames=process_fields,
        delimiter="\t",
    )

    writer.writeheader()

    for row in trace_rows:
        writer.writerow(
            {
                "sample": sample,
                "process": row.get("name", "NA"),
                "task_id": row.get("task_id", "NA"),
                "native_id": row.get("native_id", "NA"),
                "status": row.get("status", "NA"),
                "exit": row.get("exit", "NA"),
                "submit": row.get("submit", "NA"),
                "start": row.get("start", "NA"),
                "complete": row.get("complete", "NA"),
                "duration": row.get("duration", "NA"),
                "realtime": row.get("realtime", "NA"),
                "cpus": row.get("cpus", "NA"),
                "memory": row.get("memory", "NA"),
                "percent_cpu": row.get("%cpu", "NA"),
                "peak_rss": row.get("peak_rss", "NA"),
                "peak_vmem": row.get("peak_vmem", "NA"),
                "read_bytes": row.get("read_bytes", "NA"),
                "write_bytes": row.get("write_bytes", "NA"),
            }
        )


# ============================================================
# OVERALL SUMMARY
# ============================================================

completed = sum(
    1 for row in trace_rows
    if row.get("status") == "COMPLETED"
)

failed = sum(
    1 for row in trace_rows
    if row.get("status") == "FAILED"
)


summary = {
    "sample": sample,
    "process_records": len(trace_rows),
    "completed_processes": completed,
    "failed_processes": failed,
    "hardware_samples": len(hardware),

    "avg_node_cpu_temperature_C":
        f"{mean(cpu_temps):.2f}"
        if cpu_temps else "NA",

    "max_node_cpu_temperature_C":
        f"{max(cpu_temps):.2f}"
        if cpu_temps else "NA",

    "avg_node_DDR_temperature_C":
        f"{mean(ddr_temps):.2f}"
        if ddr_temps else "NA",

    "max_node_DDR_temperature_C":
        f"{max(ddr_temps):.2f}"
        if ddr_temps else "NA",

    "avg_node_RAM_used_GB":
        f"{mean(ram_values):.2f}"
        if ram_values else "NA",

    "max_node_RAM_used_GB":
        f"{max(ram_values):.2f}"
        if ram_values else "NA",

    "avg_GPU_temperature_C":
        f"{mean(gpu_temps):.2f}"
        if gpu_temps else "NA",

    "max_GPU_temperature_C":
        f"{max(gpu_temps):.2f}"
        if gpu_temps else "NA",
}


with summary_output.open("w") as fh:

    fh.write("metric\tvalue\n")

    for key, value in summary.items():
        fh.write(f"{key}\t{value}\n")


print(f"Created: {process_output}")
print(f"Created: {summary_output}")
