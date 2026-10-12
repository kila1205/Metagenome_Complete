#!/usr/bin/env bash

set -euo pipefail

echo "============================================================"
echo " METAGENOME PIPELINE - COMPUTE RESOURCE USAGE"
echo "============================================================"
echo

USER_NAME="${USER:-$(whoami)}"

echo "User: $USER_NAME"
echo "Date: $(date)"
echo

echo "============================================================"
echo " CURRENTLY RUNNING / QUEUED SLURM JOBS"
echo "============================================================"

squeue -u "$USER_NAME" \
    -o "%.12i %.10P %.35j %.8T %.10M %.6D %.6C %.12m" || true

echo
echo "============================================================"
echo " COMPLETED/RUNNING JOB RESOURCE ACCOUNTING"
echo " Last 7 days"
echo "============================================================"

sacct \
    -u "$USER_NAME" \
    -S "$(date -d '7 days ago' +%Y-%m-%d)" \
    --units=G \
    --state=RUNNING,COMPLETED,FAILED,OUT_OF_MEMORY,TIMEOUT,CANCELLED \
    --format=JobID,JobName%40,State,Elapsed,AllocCPUS,ReqMem,MaxRSS,MaxVMSize \
    -X

echo
echo "============================================================"
echo " DETAILED NEXTFLOW PROCESS RESOURCE USAGE"
echo "============================================================"

sacct \
    -u "$USER_NAME" \
    -S "$(date -d '7 days ago' +%Y-%m-%d)" \
    --units=G \
    --format=JobID,JobName%45,State,Elapsed,AllocCPUS,ReqMem,MaxRSS \
    | grep -E \
    'ASSEMB|FLYE|PRODIGAL|METABAT|GTDB|KRAKEN|BRACKEN|COVERM|MMSEQS|KOFAM|EGGNOG|REPORT|JobID' \
    || true

echo
echo "============================================================"
echo " HIGHEST MEMORY USAGE"
echo "============================================================"

sacct \
    -u "$USER_NAME" \
    -S "$(date -d '7 days ago' +%Y-%m-%d)" \
    --units=G \
    --noheader \
    --parsable2 \
    --format=JobID,JobName,State,Elapsed,AllocCPUS,ReqMem,MaxRSS \
    | awk -F'|' '
        $7 ~ /G$/ {
            x=$7
            sub(/G$/, "", x)
            if (x+0 > 0)
                print x "|" $0
        }
    ' \
    | sort -t'|' -k1,1nr \
    | head -20 \
    | cut -d'|' -f2-

echo
echo "============================================================"
echo " LONGEST RUNNING PROCESSES"
echo "============================================================"

sacct \
    -u "$USER_NAME" \
    -S "$(date -d '7 days ago' +%Y-%m-%d)" \
    --units=G \
    --noheader \
    --format=JobID,JobName%45,State,Elapsed,AllocCPUS,ReqMem,MaxRSS \
    | sort -k4 -r \
    | head -20

echo
echo "============================================================"
echo " PIPELINE RESOURCE REQUESTS FROM CONFIG"
echo "============================================================"

grep -n -E \
    'cpus|memory|time|queue|qos|clusterOptions|executor' \
    conf/cluster.config nextflow.config 2>/dev/null || true

echo
echo "============================================================"
echo " DONE"
echo "============================================================"
