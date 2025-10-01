#!/bin/bash

set -euo pipefail

PAF_FILE="../data/hprc465vschm13.aln.paf.gz"
SEQUENCE_FILES="../data/HPRC_r2_assemblies_0.6.1.agc"
REGION_PREFIX="CHM13#0#"
REFERENCE_NAME="CHM13"
THRESHOLD="0.999"
R_VALUE="5"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PICA_SCRIPT="${SCRIPT_DIR}/pica2.py"
TAJIMA_SCRIPT="${SCRIPT_DIR}/tj_d.py"

usage() {
    cat <<USAGE
Usage: $0 -b <bed_file> -l <sample_list> [options]

Required arguments:
  -b  BED file containing genomic windows
  -l  Sample list (plain text, one sequence ID per line)

Optional arguments:
  -p  PAF file for impg (default: ${PAF_FILE})
  -s  Sequence archive for impg (default: ${SEQUENCE_FILES})
  -t  Threshold for pica2.py (default: ${THRESHOLD})
  -r  R value (rounding digits) for pica2.py (default: ${R_VALUE})
  -P  Region prefix prepended to BED coordinates (default: ${REGION_PREFIX})
  -R  Reference name passed to povu gfa2vcf --stdout (default: ${REFERENCE_NAME})
  -o  Output TSV file (default: stdout)
  -h  Show this help message
USAGE
    exit 1
}

while getopts "b:l:p:s:t:r:P:R:o:h" opt; do
    case $opt in
        b) BED_FILE="$OPTARG" ;;
        l) SAMPLE_LIST="$OPTARG" ;;
        p) PAF_FILE="$OPTARG" ;;
        s) SEQUENCE_FILES="$OPTARG" ;;
        t) THRESHOLD="$OPTARG" ;;
        r) R_VALUE="$OPTARG" ;;
        P) REGION_PREFIX="$OPTARG" ;;
        R) REFERENCE_NAME="$OPTARG" ;;
        o) OUTPUT_FILE="$OPTARG" ;;
        h) usage ;;
        *) usage ;;
    esac
done

if [ -z "${BED_FILE:-}" ] || [ -z "${SAMPLE_LIST:-}" ]; then
    echo "Error: both -b <bed_file> and -l <sample_list> are required" >&2
    usage
fi

for path in "$BED_FILE" "$SAMPLE_LIST" "$PAF_FILE" "$SEQUENCE_FILES" "$PICA_SCRIPT" "$TAJIMA_SCRIPT"; do
    if [ ! -f "$path" ]; then
        echo "Error: Required file '$path' not found" >&2
        exit 1
    fi
done

for cmd in impg odgi povu python3; do
    if ! command -v "$cmd" >/dev/null 2>&1; then
        echo "Error: Required command '$cmd' not found in PATH" >&2
        exit 1
    fi
done

if ! [[ "$THRESHOLD" =~ ^[0-9]*\.?[0-9]+$ ]]; then
    echo "Error: Threshold must be numeric" >&2
    exit 1
fi

if ! [[ "$R_VALUE" =~ ^[0-9]+$ ]]; then
    echo "Error: R value must be an integer" >&2
    exit 1
fi

# Count samples (ignore blank lines and comments)
SAMPLE_COUNT=$(awk 'NF && $1 !~ /^#/' "$SAMPLE_LIST" | wc -l | awk '{print $1}')
if [ "$SAMPLE_COUNT" -lt 2 ]; then
    echo "Error: Need at least two samples to compute Tajima\'s D (found $SAMPLE_COUNT)" >&2
    exit 1
fi

if [ -n "${OUTPUT_FILE:-}" ]; then
    exec > "$OUTPUT_FILE"
fi

declare -a tmpfiles=()
cleanup() {
    for f in "${tmpfiles[@]}"; do
        [ -f "$f" ] && rm -f "$f"
    done
}
trap cleanup EXIT

printf "REGION\tLENGTH\tSAMPLES\tSEGREGATING_SITES\tPI\tTAJIMAS_D\n"

while IFS=$'\t' read -r chr start end rest; do
    if [[ -z "$chr" || "$chr" == "#"* ]]; then
        continue
    fi

    if ! [[ "$start" =~ ^[0-9]+$ && "$end" =~ ^[0-9]+$ ]]; then
        echo "Warning: Skipping malformed BED entry: $chr $start $end" >&2
        continue
    fi

    LENGTH=$((end - start))
    if [ "$LENGTH" -le 0 ]; then
        echo "Warning: Skipping non-positive interval length for $chr:$start-$end" >&2
        continue
    fi

    REGION="${REGION_PREFIX}${chr}:${start}-${end}"

    raw_gfa=$(mktemp tmp.tajd.rawgfa.XXXXXX)
    raw_og=$(mktemp tmp.tajd.og.XXXXXX)
    sorted_gfa=$(mktemp tmp.tajd.gfa.XXXXXX)
    tmpfiles+=("$raw_gfa" "$raw_og" "$sorted_gfa")

    if ! impg query -p "$PAF_FILE" -r "$REGION" --sequence-files "$SEQUENCE_FILES" -o gfa > "$raw_gfa"; then
        echo "Warning: impg query failed for region $REGION" >&2
        continue
    fi

    if ! odgi build -g "$raw_gfa" -o "$raw_og" >/dev/null 2>&1; then
        echo "Warning: odgi build failed for region $REGION" >&2
        continue
    fi

    if ! odgi sort -i "$raw_og" -o - 2>/dev/null | odgi view -i - -g > "$sorted_gfa"; then
        echo "Warning: odgi sort/view failed for region $REGION" >&2
        continue
    fi

    S_COUNT=$(povu gfa2vcf -i "$sorted_gfa" --stdout "$REFERENCE_NAME" 2>/dev/null | awk 'substr($0,1,1)!="#"' | wc -l | awk '{print $1}')
    if [ -z "$S_COUNT" ]; then
        echo "Warning: Failed to determine segregating sites for $REGION" >&2
        continue
    fi

    sim_tsv=$(mktemp tmp.tajd.sim.XXXXXX)
    tmpfiles+=("$sim_tsv")

    if ! impg similarity -p "$PAF_FILE" -r "$REGION" --sequence-files "$SEQUENCE_FILES" --subset-sequence-list "$SAMPLE_LIST" > "$sim_tsv"; then
        echo "Warning: impg similarity failed for region $REGION" >&2
        continue
    fi

    if ! pica_output=$(python3 "$PICA_SCRIPT" "$sim_tsv" -t "$THRESHOLD" -l "$LENGTH" -r "$R_VALUE" 2>/dev/null); then
        echo "Warning: pica2.py failed for region $REGION" >&2
        continue
    fi

    PI=$(printf '%s\n' "$pica_output" | awk '{print $1}')
    if [ -z "$PI" ]; then
        echo "Warning: Unable to parse pi for region $REGION" >&2
        continue
    fi

    taj_output=$(python3 "$TAJIMA_SCRIPT" -n "$SAMPLE_COUNT" -p "$PI" -S "$S_COUNT" 2>/dev/null)
    if [ -z "$taj_output" ]; then
        echo "Warning: Tajima\'s D calculation failed for region $REGION" >&2
        continue
    fi

    TAJ=$(printf '%s\n' "$taj_output" | awk '{print $3}')
    if [ -z "$TAJ" ]; then
        echo "Warning: Unable to parse Tajima\'s D for region $REGION" >&2
        continue
    fi

    if [[ "$TAJ" =~ ^nan$|^NaN$ ]]; then
        TAJ="NA"
    fi

    printf "%s\t%s\t%s\t%s\t%s\t%s\n" "$REGION" "$LENGTH" "$SAMPLE_COUNT" "$S_COUNT" "$PI" "$TAJ"

    rm -f "$raw_gfa" "$sorted_gfa" "$sim_tsv"

done < "$BED_FILE"
