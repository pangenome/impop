#!/bin/bash

set -euo pipefail

# Default paths
PAF_FILE="../data/hprc465vschm13.aln.paf.gz"
SEQUENCE_FILES="../data/HPRC_r2_assemblies_0.6.1.agc"
REGION_PREFIX="CHM13#0#"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOG_DIR="${SCRIPT_DI-o R}/fst_logs"
FST_SCRIPT="${SCRIPT_DIR}/h-fst.py"

usage() {
    cat <<USAGE
Usage: $0 -A <subset_list_A> -B <subset_list_B> -b <bed_file> [options]

Required arguments:
  -A  File with population A sequence IDs
  -B  File with population B sequence IDs
  -b  BED file containing genomic regions

Optional arguments:
  -p  PAF file for impg similarity (default: ${PAF_FILE})
  -s  Sequence files for impg similarity (default: ${SEQUENCE_FILES})
  -r  Round similarities to N decimal places (optional)
  -o  Output file (default: stdout)
  -d  Directory for log files (default: ${LOG_DIR})
  -P  Region prefix (default: ${REGION_PREFIX})
  -v  Verbose output
  -h  Display this help message

The script calculates FST = (Dxy - πxy) / Dxy for each region, where:
  - πA = nucleotide diversity within population A
  - πB = nucleotide diversity within population B
  - πxy = (πA + πB) / 2
  - Dxy = nucleotide diversity between populations
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

# Process each genomic region
process_region() {
    local chr="$1"
    local start="$2"
    local end="$3"
    local length=$((end - start))
    local region="${REGION_PREFIX}${chr}:${start}-${end}"
    
    # Create temporary file for similarities
    local tmp_sim
    tmp_sim=$(mktemp "${TMPDIR:-/tmp}/run_h-fst.similarities.XXXXXX")
    
    # Get similarities for this region
    local impg_cmd=(impg similarity -p "$PAF_FILE" -r "$region" --sequence-files "$SEQUENCE_FILES")
    
    if ! "${impg_cmd[@]}" > "$tmp_sim" 2>/dev/null; then
        echo "Error: impg similarity failed for region ${region}" >&2
        rm -f "$tmp_sim"
        return 1
    fi

    # Calculate FST
    local fst_cmd=(python3 "$FST_SCRIPT" "$tmp_sim" -a "$POP_A_FILE" -b "$POP_B_FILE" -l "$length" -d "$LOG_DIR")
    
    if [ -n "$ROUND_DIGITS" ]; then
        fst_cmd+=(-r "$ROUND_DIGITS")
    fi
    
    local fst_output
    if ! fst_output=$("${fst_cmd[@]}" 2> >(cat >&2)); then
        echo "Error: FST calculation failed for region ${region}" >&2
        rm -f "$tmp_sim"
        return 1
    fi

    # Parse output (FST, pi_A, pi_B, pi_XY, Dxy, Da)
    IFS=$'\t' read -r fst pi_a pi_b pi_xy dxy da <<< "$fst_output"
    
    # Output results
    echo -e "${region}\t${length}\t${fst}\t${pi_a}\t${pi_b}\t${pi_xy}\t${dxy}\t${da}"
    
    rm -f "$tmp_sim"
    return 0
}

# Parse command line arguments
OUTPUT_FILE=""
ROUND_DIGITS=""
VERBOSE=""

while getopts "A:B:b:p:s:r:o:d:P:vh" opt; do
    case $opt in
        A) POP_A_FILE="$OPTARG" ;;
        B) POP_B_FILE="$OPTARG" ;;
        b) BED_FILE="$OPTARG" ;;
        p) PAF_FILE="$OPTARG" ;;
        s) SEQUENCE_FILES="$OPTARG" ;;
        r) ROUND_DIGITS="$OPTARG" ;;
        o) OUTPUT_FILE="$OPTARG" ;;
        d) LOG_DIR="$OPTARG" ;;
        P) REGION_PREFIX="$OPTARG" ;;
        v) VERBOSE="1" ;;
        h) usage ;;
        *) usage ;;
    esac
done

# Validate required arguments
if [ -z "${POP_A_FILE:-}" ] || [ -z "${POP_B_FILE:-}" ] || [ -z "${BED_FILE:-}" ]; then
    echo "Error: Missing required arguments" >&2
    usage
fi

# Validate files
require_file "$POP_A_FILE" "Population A file"
require_file "$POP_B_FILE" "Population B file"
require_file "$BED_FILE" "BED file"
require_file "$PAF_FILE" "PAF file"
require_file "$SEQUENCE_FILES" "Sequence file"
require_file "$FST_SCRIPT" "fst.py script"

# Validate round digits if provided
if [ -n "$ROUND_DIGITS" ] && ! [[ "$ROUND_DIGITS" =~ ^[0-9]+$ ]]; then
    echo "Error: Round digits must be a positive integer" >&2
    exit 1
fi

# Create log directory
mkdir -p "$LOG_DIR"

# Redirect output if requested
if [ -n "$OUTPUT_FILE" ]; then
    exec > "$OUTPUT_FILE"
fi

# Output header
echo -e "REGION\tLENGTH\tFST\tPI_A\tPI_B\tPI_XY\tDXY\tDA"

# Process each region in the BED file
line_count=0
success_count=0
error_count=0

while IFS=$'\t' read -r chr start end rest; do
    ((line_count += 1))

    # Skip comment lines
    if [[ -z "$chr" || "$chr" == "#"* ]]; then
        continue
    fi

    # Validate coordinates
    if [ -z "$start" ] || [ -z "$end" ]; then
        echo "Warning: Incomplete BED entry at line ${line_count}, skipping" >&2
        ((error_count += 1))
        continue
    fi
    
    if ! [[ "$start" =~ ^[0-9]+$ && "$end" =~ ^[0-9]+$ ]]; then
        echo "Warning: Non-integer coordinates at line ${line_count}: ${chr}:${start}-${end}, skipping" >&2
        ((error_count += 1))
        continue
    fi
    
    # Validate interval
    if [ "$start" -ge "$end" ]; then
        echo "Warning: Invalid interval at line ${line_count}: ${chr}:${start}-${end}, skipping" >&2
        ((error_count += 1))
        continue
    fi
    
    # Process the region
    if [ -n "$VERBOSE" ]; then
        echo "Processing region: ${chr}:${start}-${end}" >&2
    fi
    
    if process_region "$chr" "$start" "$end"; then
        ((success_count += 1))
    else
        ((error_count += 1))
    fi
    
done < "$BED_FILE"

# Summary
if [ -n "$VERBOSE" ]; then
    echo "" >&2
    echo "Summary:" >&2
    echo "  Regions processed: ${success_count}" >&2
    echo "  Regions failed: ${error_count}" >&2
    echo "  Log directory: ${LOG_DIR}" >&2
fi
