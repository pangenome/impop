#!/bin/bash

set -euo pipefail

# Default paths
PAF_FILE="../data/hprc465vschm13.aln.paf.gz"
SEQUENCE_FILES="../data/HPRC_r2_assemblies_0.6.1.agc"
REGION_PREFIX="CHM13#0#"
LOG_DIR="./pica2_logs"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPT_PATH="${SCRIPT_DIR}/pica2.py"

usage() {
    cat <<USAGE
Usage: $0 -A <subset_list_A> -B <subset_list_B> -b <bed_file> -t <threshold> -r <r_value> [options]

Required arguments:
  -A  File with subset list A (unique sequence IDs)
  -B  File with subset list B (unique sequence IDs)
  -b  BED file containing genomic regions
  -t  Similarity threshold passed to pica2.py
  -r  R value (rounding digits) passed to pica2.py

Optional arguments:
  -p  PAF file for impg similarity (default: ${PAF_FILE})
  -s  Sequence files for impg similarity (default: ${SEQUENCE_FILES})
  -o  Write output table to file (default: stdout)
  -d  Directory to store pica2 logs (default: ${LOG_DIR})
  -P  Prefix to prepend to regions (default: ${REGION_PREFIX})
  -h  Display this help message

The script reproduces the workflow:
  1. Compute pi in subset A (piA)
  2. Compute pi in subset B (piB)
  3. Merge subset lists (A ∪ B) to form subset C and compute piC
  4. Report piA, piB, piC, average piAB = 0.5 * (piA + piB), and Fst = (piC - piAB) / piC per region
USAGE
    exit 1
}

# Input validation helper
require_file() {
    local path="$1"
    local description="$2"
    if [ ! -f "$path" ]; then
        echo "Error: ${description} '$path' not found" >&2
        exit 1
    fi
}

# Compute nucleotide diversity for a given subset
compute_pi() {
    local subset_file="$1"
    local label="$2"
    local result_var="$3"

    local tmp_sim
    tmp_sim=$(mktemp "tmp.${label}.XXXXXX")
    tmpfiles+=("$tmp_sim")

    local impg_cmd=(impg similarity -p "$PAF_FILE" -r "$REGION" --sequence-files "$SEQUENCE_FILES")
    if [ -n "$subset_file" ]; then
        impg_cmd+=(--subset-sequence-list "$subset_file")
    fi

    if ! "${impg_cmd[@]}" > "$tmp_sim"; then
        echo "Error: impg similarity failed for subset ${label} in region ${REGION}" >&2
        return 1
    fi

    local pica_output
    if ! pica_output=$(python3 "$SCRIPT_PATH" "$tmp_sim" -t "$THRESHOLD" -l "$LENGTH" -r "$R_VALUE" -d "$LOG_DIR"); then
        echo "Error: pica2.py failed for subset ${label} in region ${REGION}" >&2
        echo "$pica_output" >&2
        return 1
    fi

    local pi_value
    pi_value=$(printf '%s\n' "$pica_output" | awk '{print $1}')
    if [ -z "$pi_value" ]; then
        echo "Error: Unable to parse pi value for subset ${label} in region ${REGION}" >&2
        return 1
    fi

    printf -v "$result_var" '%s' "$pi_value"

    rm -f "$tmp_sim"
    return 0
}

OUTPUT_FILE=""

while getopts "A:B:b:p:s:t:r:o:d:P:h" opt; do
    case $opt in
        A) SUBSET_A_LIST="$OPTARG" ;;
        B) SUBSET_B_LIST="$OPTARG" ;;
        b) BED_FILE="$OPTARG" ;;
        p) PAF_FILE="$OPTARG" ;;
        s) SEQUENCE_FILES="$OPTARG" ;;
        t) THRESHOLD="$OPTARG" ;;
        r) R_VALUE="$OPTARG" ;;
        o) OUTPUT_FILE="$OPTARG" ;;
        d) LOG_DIR="$OPTARG" ;;
        P) REGION_PREFIX="$OPTARG" ;;
        h) usage ;;
        *) usage ;;
    esac
done

# Ensure required arguments are present
if [ -z "${SUBSET_A_LIST:-}" ] || [ -z "${SUBSET_B_LIST:-}" ] || [ -z "${BED_FILE:-}" ] || [ -z "${THRESHOLD:-}" ] || [ -z "${R_VALUE:-}" ]; then
    echo "Error: Missing required arguments" >&2
    usage
fi

# Validate files
require_file "$SUBSET_A_LIST" "Subset list A"
require_file "$SUBSET_B_LIST" "Subset list B"
require_file "$BED_FILE" "BED file"
require_file "$PAF_FILE" "PAF file"
require_file "$SEQUENCE_FILES" "Sequence file"
require_file "$SCRIPT_PATH" "pica2.py script"

# Validate threshold and R value
if ! [[ "$THRESHOLD" =~ ^[0-9]*\.?[0-9]+$ ]]; then
    echo "Error: Threshold must be numeric" >&2
    exit 1
fi

if ! [[ "$R_VALUE" =~ ^[0-9]+$ ]]; then
    echo "Error: R value must be an integer" >&2
    exit 1
fi

# Prepare log directory
mkdir -p "$LOG_DIR"

# Track temporary files for cleanup
trap 'for f in "${tmpfiles[@]:-}"; do [ -f "$f" ] && rm -f "$f"; done' EXIT

declare -a tmpfiles=()

# Build union list (subset C)
UNION_LIST=$(mktemp tmp.union.XXXXXX)
tmpfiles+=("$UNION_LIST")
awk 'NF && !seen[$0]++' "$SUBSET_A_LIST" "$SUBSET_B_LIST" > "$UNION_LIST"

if [ ! -s "$UNION_LIST" ]; then
    echo "Error: Union subset list is empty" >&2
    exit 1
fi

SUBSET_A_LABEL=$(basename "$SUBSET_A_LIST")
SUBSET_B_LABEL=$(basename "$SUBSET_B_LIST")
SUBSET_C_LABEL=$(basename "$SUBSET_A_LIST")"+"$(basename "$SUBSET_B_LIST")

# Redirect output if requested
if [ -n "$OUTPUT_FILE" ]; then
    exec > "$OUTPUT_FILE"
fi

echo -e "REGION\tLENGTH\tTHRESHOLD\tR_VALUE\tPI_A\tPI_B\tPI_C\tPI_AB_AVG\tFST"

# Process each region in the BED file
while IFS=$'\t' read -r chr start end rest; do
    if [[ -z "$chr" || "$chr" == "#"* ]]; then
        continue
    fi

    if [ -z "$start" ] || [ -z "$end" ]; then
        echo "Warning: Incomplete BED entry for chromosome ${chr}, skipping" >&2
        continue
    fi

    if ! [[ "$start" =~ ^[0-9]+$ && "$end" =~ ^[0-9]+$ ]]; then
        echo "Warning: Non-integer coordinates in BED entry ${chr}\t${start}\t${end}, skipping" >&2
        continue
    fi

    LENGTH=$((end - start))
    if [ "$LENGTH" -le 0 ]; then
        echo "Warning: Non-positive interval length for ${chr}:${start}-${end}, skipping" >&2
        continue
    fi

    REGION="${REGION_PREFIX}${chr}:${start}-${end}"

    if ! compute_pi "$SUBSET_A_LIST" "A" PI_A; then
        echo "Warning: Skipping region ${REGION} due to errors with subset A" >&2
        continue
    fi

    if ! compute_pi "$SUBSET_B_LIST" "B" PI_B; then
        echo "Warning: Skipping region ${REGION} due to errors with subset B" >&2
        continue
    fi

    if ! compute_pi "$UNION_LIST" "C" PI_C; then
        echo "Warning: Skipping region ${REGION} due to errors with subset C" >&2
        continue
    fi

    PI_AB_AVG=$(python3 - "$PI_A" "$PI_B" <<'PY'
import sys
pi_a = float(sys.argv[1])
pi_b = float(sys.argv[2])
print(f"{0.5 * (pi_a + pi_b):.8f}")
PY
    )

    FST=$(python3 - "$PI_A" "$PI_B" "$PI_C" <<'PY'
import sys
pi_a = float(sys.argv[1])
pi_b = float(sys.argv[2])
pi_c = float(sys.argv[3])
pi_ab = 0.5 * (pi_a + pi_b)
if pi_c == 0:
    print("NA")
else:
    print(f"{(pi_c - pi_ab) / pi_c:.8f}")
PY
    )

    echo -e "${REGION}\t${LENGTH}\t${THRESHOLD}\t${R_VALUE}\t${PI_A}\t${PI_B}\t${PI_C}\t${PI_AB_AVG}\t${FST}"

done < "$BED_FILE"
