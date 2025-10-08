#!/usr/bin/env bash

# Run Tajima's D calculations for several population panels.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
RUN_TAJD="${PROJECT_ROOT}/scripts/run_tajd.sh"
PAF_FILE="${PROJECT_ROOT}/data/hprc465vschm13.aln.paf.gz"
SEQUENCE_FILES="${PROJECT_ROOT}/data/HPRC_r2_assemblies_0.6.1.agc"
METADATA_DIR="${PROJECT_ROOT}/../metadata"

usage() {
    cat <<USAGE
Usage: $(basename "$0") [-b bed_file] [-p paf_file] [-s agc_file]

Options:
  -b  BED file listing regions to process (default: region.bed in current directory)
  -p  Override the default PAF file (default: ${PAF_FILE})
  -s  Override the default AGC file (default: ${SEQUENCE_FILES})
  -h  Show this help text and exit

Outputs are written alongside your current working directory.
USAGE
}

BED_FILE="region.bed"

while getopts ":b:p:s:h" opt; do
    case "$opt" in
        b) BED_FILE="$OPTARG" ;;
        p) PAF_FILE="$OPTARG" ;;
        s) SEQUENCE_FILES="$OPTARG" ;;
        h) usage; exit 0 ;;
        :) echo "Error: -$OPTARG requires a value" >&2; usage; exit 1 ;;
        *) usage; exit 1 ;;
    esac
done

if [ ! -f "$RUN_TAJD" ]; then
    echo "Error: run_tajd.sh not found at $RUN_TAJD" >&2
    exit 1
fi

if [ ! -f "$PAF_FILE" ]; then
    echo "Error: PAF file not found at $PAF_FILE" >&2
    exit 1
fi

if [ ! -f "$SEQUENCE_FILES" ]; then
    echo "Error: sequence file not found at $SEQUENCE_FILES" >&2
    exit 1
fi

if [ ! -f "$BED_FILE" ]; then
    echo "Error: BED file not found at $BED_FILE" >&2
    exit 1
fi

panels=(
    "EUR eur.tj"
    "AFR afr.tj"
    "EAS eas.tj"
    "SAS sas.tj"
    "AMR amr.tj"
)

for entry in "${panels[@]}"; do
    read -r group output_file <<<"$entry"
    subset_file="${METADATA_DIR}/agc.${group}"

    if [ ! -f "$subset_file" ]; then
        echo "Error: subset list not found: $subset_file" >&2
        exit 1
    fi

    echo "[tajd] ${group} -> ${output_file}"
    "${RUN_TAJD}" \
        -b "$BED_FILE" \
        -l "$subset_file" \
        -p "$PAF_FILE" \
        -s "$SEQUENCE_FILES" \
        -o "$output_file"
done
