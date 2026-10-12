#!/usr/bin/env bash

set -u

OUTPUT_FILE="${1:?Usage: monitor_hardware.sh OUTPUT_FILE [INTERVAL]}"
INTERVAL="${2:-60}"

mkdir -p "$(dirname "$OUTPUT_FILE")"

NODE="$(hostname)"

printf '%s\n' \
"timestamp	epoch	node	load_1m	load_5m	load_15m	cpu_mhz_avg	cpu_mhz_max	ram_total_GB	ram_used_GB	ram_available_GB	swap_used_GB	CPU1_temp_C	CPU2_temp_C	DDR4_A_temp_C	DDR4_B_temp_C	DDR4_C_temp_C	DDR4_D_temp_C	DDR4_I_temp_C	DDR4_J_temp_C	DDR4_K_temp_C	DDR4_L_temp_C	VOLT_12V	VOLT_CPU1_VCORE	VOLT_CPU2_VCORE	FAN1_RPM	FAN2_RPM	FAN3_RPM	FAN4_RPM	GPU0_temp_C	GPU0_power_W	GPU0_util_pct	GPU0_mem_used_MB	GPU1_temp_C	GPU1_power_W	GPU1_util_pct	GPU1_mem_used_MB	GPU2_temp_C	GPU2_power_W	GPU2_util_pct	GPU2_mem_used_MB" \
> "$OUTPUT_FILE"


get_ipmi_value() {
    local sensor="$1"
    local data="$2"
    local value

    value="$(
        printf '%s\n' "$data" |
        awk -F'|' -v sensor="$sensor" '
        $1 ~ sensor {
            gsub(/^[ \t]+|[ \t]+$/, "", $2)
            split($2,a," ")
            print a[1]
            exit
        }'
    )"

    if [[ -z "$value" ]]; then
        echo "NA"
    else
        echo "$value"
    fi
}


get_gpu() {
    local gpu_data="$1"
    local line_number="$2"
    local line

    line="$(printf '%s\n' "$gpu_data" | sed -n "${line_number}p")"

    if [[ -z "$line" ]]; then
        printf 'NA\tNA\tNA\tNA'
        return
    fi

    printf '%s\n' "$line" |
    awk -F',' '
    {
        for(i=1;i<=4;i++) {
            gsub(/^[ \t]+|[ \t]+$/, "", $i)
            if($i == "" || $i == "[N/A]") {
                $i="NA"
            }
        }
        printf "%s\t%s\t%s\t%s",$1,$2,$3,$4
    }'
}


cleanup() {
    exit 0
}

trap cleanup SIGTERM SIGINT


while true
do
    TIMESTAMP="$(date '+%Y-%m-%d %H:%M:%S')"
    EPOCH="$(date +%s)"

    read -r LOAD1 LOAD5 LOAD15 _ < /proc/loadavg

    CPU_AVG="$(
        awk '
        /cpu MHz/ {
            sum += $4
            n++
        }
        END {
            if(n>0) printf "%.0f",sum/n
            else print "NA"
        }' /proc/cpuinfo
    )"

    CPU_MAX="$(
        awk '
        /cpu MHz/ {
            if($4>max) max=$4
        }
        END {
            if(max>0) printf "%.0f",max
            else print "NA"
        }' /proc/cpuinfo
    )"

    read -r RAM_TOTAL RAM_USED RAM_AVAILABLE SWAP_USED < <(
        free -m |
        awk '
        /^Mem:/ {
            total=$2/1024
            used=$3/1024
            available=$7/1024
        }
        /^Swap:/ {
            swap=$3/1024
        }
        END {
            printf "%.2f %.2f %.2f %.2f\n",
                total,used,available,swap
        }'
    )

    IPMI=""

    if command -v ipmitool >/dev/null 2>&1; then
        IPMI="$(ipmitool sensor 2>/dev/null || true)"
    fi

    CPU1="$(get_ipmi_value "TEMP_CPU1" "$IPMI")"
    CPU2="$(get_ipmi_value "TEMP_CPU2" "$IPMI")"

    DDRA="$(get_ipmi_value "TEMP_DDR4_A" "$IPMI")"
    DDRB="$(get_ipmi_value "TEMP_DDR4_B" "$IPMI")"
    DDRC="$(get_ipmi_value "TEMP_DDR4_C" "$IPMI")"
    DDRD="$(get_ipmi_value "TEMP_DDR4_D" "$IPMI")"

    DDRI="$(get_ipmi_value "TEMP_DDR4_I" "$IPMI")"
    DDRJ="$(get_ipmi_value "TEMP_DDR4_J" "$IPMI")"
    DDRK="$(get_ipmi_value "TEMP_DDR4_K" "$IPMI")"
    DDRL="$(get_ipmi_value "TEMP_DDR4_L" "$IPMI")"

    V12="$(get_ipmi_value "VOLT_12V" "$IPMI")"
    VCORE1="$(get_ipmi_value "VOLT_CPU1_VCORE" "$IPMI")"
    VCORE2="$(get_ipmi_value "VOLT_CPU2_VCORE" "$IPMI")"

    FAN1="$(get_ipmi_value "FAN1_1" "$IPMI")"
    FAN2="$(get_ipmi_value "FAN2_1" "$IPMI")"
    FAN3="$(get_ipmi_value "FAN3_1" "$IPMI")"
    FAN4="$(get_ipmi_value "FAN4_1" "$IPMI")"

    GPU_DATA=""

    if command -v nvidia-smi >/dev/null 2>&1; then
        GPU_DATA="$(
            nvidia-smi \
            --query-gpu=temperature.gpu,power.draw,utilization.gpu,memory.used \
            --format=csv,noheader,nounits \
            2>/dev/null || true
        )"
    fi

    GPU0="$(get_gpu "$GPU_DATA" 1)"
    GPU1="$(get_gpu "$GPU_DATA" 2)"
    GPU2="$(get_gpu "$GPU_DATA" 3)"

    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$TIMESTAMP" \
        "$EPOCH" \
        "$NODE" \
        "$LOAD1" \
        "$LOAD5" \
        "$LOAD15" \
        "$CPU_AVG" \
        "$CPU_MAX" \
        "$RAM_TOTAL" \
        "$RAM_USED" \
        "$RAM_AVAILABLE" \
        "$SWAP_USED" \
        "$CPU1" \
        "$CPU2" \
        "$DDRA" \
        "$DDRB" \
        "$DDRC" \
        "$DDRD" \
        "$DDRI" \
        "$DDRJ" \
        "$DDRK" \
        "$DDRL" \
        "$V12" \
        "$VCORE1" \
        "$VCORE2" \
        "$FAN1" \
        "$FAN2" \
        "$FAN3" \
        "$FAN4" \
        $GPU0 \
        $GPU1 \
        $GPU2 \
        >> "$OUTPUT_FILE"

    sleep "$INTERVAL"
done
