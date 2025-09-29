#!/bin/bash

set -euo pipefail

# Default paths
PAF_FILE="../data/hprc465vschm13.aln.paf.gz"
SEQUENCE_FILES="../data/HPRC_r2_assemblies_0.6.1.agc"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPT_PATH="${SCRIPT_DIR}/pica2.py"

# Display usage information
usage() {
    cat <<USAGE
Usage: $0 -b <bed_file> -t <threshold> -r <r_value> [options]
  -b  BED file containing regions (required)
  -p  PAF file for impg similarity (default: ${PAF_FILE})
  -s  Sequence files for impg similarity (default: ${SEQUENCE_FILES})
  -u  File with assemblies to subset (passed to --subset-sequence-list)

  pica2 options:
    -t  Similarity threshold for pica2.py (required)
    -r  R value for pica2.py (required)

  -h  Display this help message

Example:
  $0 -b ackr1.win.bed -t 0.999 -r 4 -u ../metadata/agc.EUR
USAGE
    exit 1
}

# Parse command line options
while getopts "b:t:r:p:s:u:h" opt; do
    case $opt in
        b) BED_FILE="$OPTARG" ;;
        t) THRESHOLD="$OPTARG" ;;
        r) R_VALUE="$OPTARG" ;;
        p) PAF_FILE="$OPTARG" ;;
        s) SEQUENCE_FILES="$OPTARG" ;;
        u) SUBSET_LIST="$OPTARG" ;;
        h) usage ;;
        *) usage ;;
    esac
done

# Ensure required arguments are present
if [ -z "${BED_FILE:-}" ] || [ -z "${THRESHOLD:-}" ] || [ -z "${R_VALUE:-}" ]; then
    echo "Error: Missing required arguments" >&2
    usage
fi

# Validate files
if [ ! -f "$BED_FILE" ]; then
    echo "Error: BED file '$BED_FILE' not found" >&2
    exit 1
fi

if [ ! -f "$PAF_FILE" ]; then
    echo "Error: PAF file '$PAF_FILE' not found" >&2
    exit 1
fi

if [ ! -f "$SEQUENCE_FILES" ]; then
    echo "Error: Sequence file '$SEQUENCE_FILES' not found" >&2
    exit 1
fi

if [ -n "${SUBSET_LIST:-}" ] && [ ! -f "$SUBSET_LIST" ]; then
    echo "Error: Subset list '$SUBSET_LIST' not found" >&2
    exit 1
fi

# Validate threshold (float)
if ! [[ "$THRESHOLD" =~ ^[0-9]*\.?[0-9]+$ ]]; then
    echo "Error: Threshold must be a number (e.g., 0.988)" >&2
    exit 1
fi

# Validate R value (integer)
if ! [[ "$R_VALUE" =~ ^[0-9]+$ ]]; then
    echo "Error: R value must be an integer (e.g., 5)" >&2
    exit 1
fi

# Ensure temporary files are removed even on failure
tmpfiles=()
cleanup() {
    for f in "${tmpfiles[@]}"; do
        [ -f "$f" ] && rm -f "$f"
    done
}
trap cleanup EXIT

# Print table header
echo -e "REGION\tLENGTH\tTHRESHOLD\tR_VALUE\tPICA_OUTPUT"

# Process each window in the BED file
while IFS=$'\t' read -r chr start end name; do
    # Skip comments or empty rows
    if [[ -z "$chr" || "$chr" == "#"* ]]; then
        continue
    fi

    # Compute window length
    LENGTH=$((end - start))
    if [ "$LENGTH" -le 0 ]; then
        echo "Warning: Skipping region with non-positive length: $chr:$start-$end" >&2
        continue
    fi

    REGION="CHM13#0#${chr}:${start}-${end}"
    tmp_sim=$(mktemp tmp.sim.XXXXXX)
    tmpfiles+=("$tmp_sim")

    # Build impg command
    impg_cmd=(impg similarity -p "$PAF_FILE" -r "$REGION" --sequence-files "$SEQUENCE_FILES")
    if [ -n "${SUBSET_LIST:-}" ]; then
        impg_cmd+=(--subset-sequence-list "$SUBSET_LIST")
    fi

    # Step 1: Generate similarity matrix via impg
    if ! "${impg_cmd[@]}" > "$tmp_sim" 2>/dev/null; then
        echo "Error: impg similarity failed for region $REGION" >&2
        continue
    fi

    # Step 2: Evaluate nucleotide diversity for the window
    if ! PICA_RAW=$(python3 "$SCRIPT_PATH" "$tmp_sim" -t "$THRESHOLD" -l "$LENGTH" -r "$R_VALUE" 2>&1); then
        echo "Error: pica2.py failed for region $REGION" >&2
        echo "$PICA_RAW" >&2
        continue
    fi

    PICA_OUTPUT=$(echo "$PICA_RAW" | tr '\n' ' ' | sed 's/  */ /g' | sed 's/ *$//')

    echo -e "${REGION}\t${LENGTH}\t${THRESHOLD}\t${R_VALUE}\t${PICA_OUTPUT}"

done < "$BED_FILE"
