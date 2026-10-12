#!/usr/bin/env bash

set -euo pipefail

TRACE="${1:?Usage: collect_slurm_resources.sh TRACE_FILE OUTPUT_FILE}"
OUTPUT="${2:?Usage: collect_slurm_resources.sh TRACE_FILE OUTPUT_FILE}"

if [[ ! -s "$TRACE" ]]
then
    echo "ERROR: Trace file does not exist or is empty: $TRACE" >&2
    exit 1
fi

mkdir -p "$(dirname "$OUTPUT")"

# ============================================================
# LOCATE native_id COLUMN
# ============================================================

HEADER="$(head -1 "$TRACE")"

NATIVE_COL="$(
    printf '%s\n' "$HEADER" |
    awk -F'\t' '
    {
        for(i=1;i<=NF;i++) {
            if($i=="native_id") {
                print i
                exit
            }
        }
    }'
)"

if [[ -z "$NATIVE_COL" ]]
then
    echo "ERROR: native_id column not found in trace." >&2
    exit 1
fi

# ============================================================
# EXTRACT UNIQUE SLURM JOB IDs
# ============================================================

mapfile -t JOB_IDS < <(
    tail -n +2 "$TRACE" |
    awk -F'\t' -v col="$NATIVE_COL" '
        $col != "" && $col != "-" {
            print $col
        }
    ' |
    sort -u
)

if [[ "${#JOB_IDS[@]}" -eq 0 ]]
then
    echo "ERROR: No Slurm native IDs found in trace." >&2
    exit 1
fi

JOB_LIST="$(
    printf '%s\n' "${JOB_IDS[@]}" |
    paste -sd, -
)"

echo "Found ${#JOB_IDS[@]} Slurm jobs."
echo "Job IDs: $JOB_LIST"

# ============================================================
# QUERY EXACT SLURM JOBS
# ============================================================

TMP="$(mktemp)"

sacct \
    -j "$JOB_LIST" \
    --units=G \
    --parsable2 \
    --noheader \
    --format=JobID,JobName,State,ExitCode,NodeList,AllocCPUS,ReqMem,Elapsed,TotalCPU,CPUTime,MaxRSS,MaxVMSize,AveRSS,AveCPU,MaxDiskRead,MaxDiskWrite \
    > "$TMP"

# sacct --parsable2 is pipe-delimited.
# Convert it into a true TSV.

{
    printf 'JobID\tJobName\tState\tExitCode\tNodeList\tAllocCPUS\tReqMem\tElapsed\tTotalCPU\tCPUTime\tMaxRSS\tMaxVMSize\tAveRSS\tAveCPU\tMaxDiskRead\tMaxDiskWrite\n'

    tr '|' '\t' < "$TMP"

} > "$OUTPUT"

rm -f "$TMP"

echo "Created: $OUTPUT"
