#!/usr/bin/env bash

# Run a batch of h-fst comparisons across population panels.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
RUN_H_FST="${PROJECT_ROOT}/scripts/run_h-fst.sh"
PAF_FILE="${PROJECT_ROOT}/../data/hprc465vschm13.aln.paf.gz"
SEQUENCE_FILES="${PROJECT_ROOT}/../data/HPRC_r2_assemblies_0.6.1.agc"
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

if [ ! -f "$RUN_H_FST" ]; then
    echo "Error: run_h-fst.sh not found at $RUN_H_FST" >&2
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

pairs=(
    "EUR AFR eur.afr.fst"
    "EAS AFR eas.afr.fst"
    "SAS AFR sas.afr.fst"
    "AMR AFR amr.afr.fst"
    "EAS EUR eas.eur.fst"
    "SAS EUR sas.eur.fst"
    "AMR EUR amr.eur.fst"
    "EAS SAS eas.sas.fst"
    "AMR SAS amr.sas.fst"
    "AMR EAS amr.eas.fst"
)

for entry in "${pairs[@]}"; do
    read -r group_a group_b output_file <<<"$entry"
    subset_a="${METADATA_DIR}/agc.${group_a}"
    subset_b="${METADATA_DIR}/agc.${group_b}"

    if [ ! -f "$subset_a" ]; then
        echo "Error: subset list not found: $subset_a" >&2
        exit 1
    fi
    if [ ! -f "$subset_b" ]; then
        echo "Error: subset list not found: $subset_b" >&2
        exit 1
    fi

    echo "[h-fst] ${group_a} vs ${group_b} -> ${output_file}"
    "${RUN_H_FST}" \
        -p "$PAF_FILE" \
        -s "$SEQUENCE_FILES" \
        -b "$BED_FILE" \
        -A "$subset_a" \
        -B "$subset_b" \
        -o "$output_file"
done
