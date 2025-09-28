#!/bin/bash

set -euo pipefail

# Default paths
PAF_FILE="../data/hprc465vschm13.aln.paf.gz"
SEQUENCE_FILES="../data/HPRC_r2_assemblies_0.6.1.agc"
SCRIPT_PATH="pica2.py"

# Display usage information
usage() {
    echo "Usage: $0 -b <bed_file> -t <threshold> -r <r_value>"
    echo "  -b  BED file containing regions (required)"
    echo "  -t  Similarity threshold for pica2.2.py (required)"
    echo "  -r  R value for pica2.2.py (required)"
    echo "  -h  Display this help message"
    echo ""
    echo "Example: $0 -b ackr1_win.bed -t 0.988 -r 5"
    exit 1
}

# Parse command line options
while getopts "b:t:r:h" opt; do
    case $opt in
        b) BED_FILE="$OPTARG" ;;
        t) THRESHOLD="$OPTARG" ;;
        r) R_VALUE="$OPTARG" ;;
        h) usage ;;
        *) usage ;;
    esac
done

# Ensure required arguments are present
if [ -z "${BED_FILE:-}" ] || [ -z "${THRESHOLD:-}" ] || [ -z "${R_VALUE:-}" ]; then
    echo "Error: Missing required arguments"
    usage
fi

# Validate BED file exists
if [ ! -f "$BED_FILE" ]; then
    echo "Error: BED file '$BED_FILE' not found"
    exit 1
fi

# Validate threshold (float)
if ! [[ "$THRESHOLD" =~ ^[0-9]*\.?[0-9]+$ ]]; then
    echo "Error: Threshold must be a number (e.g., 0.988)"
    exit 1
fi

# Validate R value (integer)
if ! [[ "$R_VALUE" =~ ^[0-9]+$ ]]; then
    echo "Error: R value must be an integer (e.g., 5)"
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

    # Step 1: Generate similarity matrix directly via impg
    if ! impg similarity -p "$PAF_FILE" -r "$REGION" --sequence-files "$SEQUENCE_FILES" > "$tmp_sim" 2>/dev/null; then
        echo "Error: impg similarity failed for region $REGION" >&2
        continue
    fi

    # Step 2: Evaluate nucleotide diversity for the window
    PICA_OUTPUT=$(python3 "$SCRIPT_PATH" "$tmp_sim" -t "$THRESHOLD" -l "$LENGTH" -r "$R_VALUE" 2>&1 | tr '\n' ' ' | sed 's/  */ /g' | sed 's/ *$//')

    echo -e "${REGION}\t${LENGTH}\t${THRESHOLD}\t${R_VALUE}\t${PICA_OUTPUT}"

done < "$BED_FILE"
