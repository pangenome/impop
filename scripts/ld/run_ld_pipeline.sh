#!/usr/bin/env bash
# Orchestrate multi-locus LD analysis for all region pairs provided in a BED file.

set -euo pipefail

usage() {
    cat <<'USAGE'
Usage:
  run_ld_pipeline.sh -b regions.bed -p alignments.paf -f sequences.fasta -n SAMPLE_SIZE [options]
  run_ld_pipeline.sh --region CHR:START-END --interval-size BP -p alignments.paf -f sequences.fasta -n SAMPLE_SIZE [options]

Required arguments:
  -p, --paf FILE            Path to the wfmash PAF file
  -f, --fasta FILE          Reference FASTA passed to impg
  -n, --sample-size INT     Sample size N used for LD summaries

Region input (choose one mode):
  -b, --bed FILE            BED file with >=2 loci (optional 4th column for locus name)
      --region RANGE        Genomic range in the format CHR:START-END
      --interval-size BP    Interval size in bp used to tile --region

Optional arguments:
  -o, --out-dir DIR         Directory for outputs (default: current directory)
      --impg PATH           Path to the impg executable (default: impg)
      --subset-sequence-list FILE
                             Restrict impg similarity to sequence names listed in FILE
      --max-distance BP     Maximum allowed locus separation in bp (default: 1000000)
      --progress            Print progress messages to stderr
      --debug               Keep pair-level intermediate files and directories
  -h, --help                Show this help message and exit

Outputs:
  * OUT_DIR/ld_summary.tsv with one row per evaluated pair
  * OUT_DIR/ld_heatmap.tsv lower-triangular matrix (Cramer's V)
  * OUT_DIR/loci/*.sim cached once per locus
  * With --debug, per-pair intermediate files under OUT_DIR/<locusA>_<locusB>/
USAGE
}

BED_FILE=""
REGION_RANGE=""
INTERVAL_SIZE=""
PAF_FILE=""
FASTA_FILE=""
SAMPLE_SIZE=""
OUT_DIR="."
IMPG_BIN="impg"
SUBSET_SEQUENCE_LIST=""
MAX_DISTANCE=1000000
DEBUG_MODE=0
SHOW_PROGRESS=0

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
CLUSTER_SCRIPT="$SCRIPT_DIR/cluster_two_loci.py"
DIJ_SCRIPT="$SCRIPT_DIR/compute_dij.py"

log() {
    if (( SHOW_PROGRESS )); then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" >&2
    fi
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        -b|--bed)
            BED_FILE=$2
            shift 2
            ;;
        --region)
            REGION_RANGE=$2
            shift 2
            ;;
        --interval-size)
            INTERVAL_SIZE=$2
            shift 2
            ;;
        -p|--paf)
            PAF_FILE=$2
            shift 2
            ;;
        -f|--fasta|--sequence-files)
            FASTA_FILE=$2
            shift 2
            ;;
        -n|--sample-size)
            SAMPLE_SIZE=$2
            shift 2
            ;;
        -o|--out-dir)
            OUT_DIR=$2
            shift 2
            ;;
        --impg)
            IMPG_BIN=$2
            shift 2
            ;;
        --subset-sequence-list)
            SUBSET_SEQUENCE_LIST=$2
            shift 2
            ;;
        --max-distance)
            MAX_DISTANCE=$2
            shift 2
            ;;
        --debug)
            DEBUG_MODE=1
            shift
            ;;
        --progress)
            SHOW_PROGRESS=1
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "Unknown argument: $1" >&2
            usage
            exit 1
            ;;
    esac
done

if [[ -z $PAF_FILE || -z $FASTA_FILE || -z $SAMPLE_SIZE ]]; then
    echo "Error: --paf, --fasta, and --sample-size are required." >&2
    usage
    exit 1
fi

if [[ -n $BED_FILE && -n $REGION_RANGE ]]; then
    echo "Error: provide either --bed or --region/--interval-size, not both." >&2
    usage
    exit 1
fi

if [[ -z $BED_FILE && -z $REGION_RANGE ]]; then
    echo "Error: either --bed or --region must be provided." >&2
    usage
    exit 1
fi

if [[ -n $BED_FILE && -n $INTERVAL_SIZE ]]; then
    echo "Error: --interval-size can only be used with --region." >&2
    usage
    exit 1
fi

if [[ -n $REGION_RANGE && -z $INTERVAL_SIZE ]]; then
    echo "Error: --interval-size is required when using --region." >&2
    usage
    exit 1
fi

if [[ ! -f $CLUSTER_SCRIPT ]]; then
    echo "Error: cluster_two_loci.py not found next to this script." >&2
    exit 1
fi

if [[ ! -f $DIJ_SCRIPT ]]; then
    echo "Error: compute_dij.py not found next to this script." >&2
    exit 1
fi

if ! command -v "$IMPG_BIN" >/dev/null 2>&1; then
    echo "Error: impg executable not found (looked for '$IMPG_BIN')." >&2
    exit 1
fi

if [[ -n $BED_FILE && ! -s $BED_FILE ]]; then
    echo "Error: BED file '$BED_FILE' not found or empty." >&2
    exit 1
fi

if [[ -n $SUBSET_SEQUENCE_LIST && ! -s $SUBSET_SEQUENCE_LIST ]]; then
    echo "Error: subset sequence list '$SUBSET_SEQUENCE_LIST' not found or empty." >&2
    exit 1
fi

if ! [[ $MAX_DISTANCE =~ ^[0-9]+$ ]]; then
    echo "Error: --max-distance must be an integer number of base pairs." >&2
    exit 1
fi

if ! [[ $SAMPLE_SIZE =~ ^[0-9]+(\.[0-9]+)?$ ]]; then
    echo "Error: --sample-size must be numeric." >&2
    exit 1
fi

if [[ -n $INTERVAL_SIZE && ! $INTERVAL_SIZE =~ ^[0-9]+$ ]]; then
    echo "Error: --interval-size must be an integer number of base pairs." >&2
    exit 1
fi

declare -a REFS STARTS ENDS LABELS
idx=0
if [[ -n $BED_FILE ]]; then
    while IFS=$' \t' read -r ref start end label rest; do
        [[ -z $ref ]] && continue
        [[ $ref =~ ^# ]] && continue
        if [[ -z $start || -z $end ]]; then
            log "Skipping malformed BED line: $ref $start $end"
            continue
        fi
        if [[ -z $label ]]; then
            label="loc$((idx + 1))"
        fi
        REFS[idx]=$ref
        STARTS[idx]=$start
        ENDS[idx]=$end
        LABELS[idx]=$label
        idx=$((idx + 1))
    done < "$BED_FILE"
else
    if [[ ! $REGION_RANGE =~ ^([^:]+):([0-9]+)-([0-9]+)$ ]]; then
        echo "Error: --region must match CHR:START-END." >&2
        exit 1
    fi

    region_ref=${BASH_REMATCH[1]}
    region_start=${BASH_REMATCH[2]}
    region_end=${BASH_REMATCH[3]}

    if (( region_end <= region_start )); then
        echo "Error: --region end must be greater than start." >&2
        exit 1
    fi

    current_start=$region_start
    while (( current_start < region_end )); do
        current_end=$((current_start + INTERVAL_SIZE))
        if (( current_end > region_end )); then
            current_end=$region_end
        fi
        REFS[idx]=$region_ref
        STARTS[idx]=$current_start
        ENDS[idx]=$current_end
        LABELS[idx]="loc$((idx + 1))"
        idx=$((idx + 1))
        current_start=$current_end
    done
fi

NUM_LOCI=${#REFS[@]}
if [[ $NUM_LOCI -lt 2 ]]; then
    echo "Error: fewer than two valid regions found in BED file." >&2
    exit 1
fi

log "Loaded $NUM_LOCI loci from $BED_FILE"

if (( NUM_LOCI > 100 )); then
    log "Warning: BED file contains $NUM_LOCI loci; pairwise LD will evaluate many combinations."
fi

declare -A MIN_STARTS MAX_ENDS
for ((i=0; i<NUM_LOCI; i++)); do
    ref=${REFS[i]}
    start=${STARTS[i]}
    end=${ENDS[i]}

    if [[ -z ${MIN_STARTS[$ref]+x} || start -lt MIN_STARTS[$ref] ]]; then
        MIN_STARTS[$ref]=$start
    fi
    if [[ -z ${MAX_ENDS[$ref]+x} || end -gt MAX_ENDS[$ref] ]]; then
        MAX_ENDS[$ref]=$end
    fi
done

for ref in "${!MIN_STARTS[@]}"; do
    span=$((MAX_ENDS[$ref] - MIN_STARTS[$ref]))
    if (( span > 5000000 )); then
        log "Warning: BED loci span $span bp on $ref; disk usage can grow substantially, especially with --debug."
    fi
done

mkdir -p "$OUT_DIR"
SUMMARY_PATH="$OUT_DIR/ld_summary.tsv"
HEATMAP_PATH="$OUT_DIR/ld_heatmap.tsv"
LOCI_DIR="$OUT_DIR/loci"
TMP_STATS=$(mktemp)
trap 'rm -f "$TMP_STATS"' EXIT

printf "ref\tl1start\tl1end\tl2ref\tl2start\tl2end\tsample_size\tchi_squared\tcramers_v\tcomposite_d_prime\tmultiallelic_r2\n" > "$SUMMARY_PATH"

declare -A HEATMAP
declare -a SIM_PATHS
for ((i=0; i<NUM_LOCI; i++)); do
    HEATMAP["$i,$i"]=1.0
done

mkdir -p "$LOCI_DIR"

for ((i=0; i<NUM_LOCI; i++)); do
    ref=${REFS[i]}
    start=${STARTS[i]}
    end=${ENDS[i]}
    label=${LABELS[i]}

    sim_path="$LOCI_DIR/${label}.sim"
    SIM_PATHS[i]=$sim_path

    log "[$label] similarity for $ref:$start-$end"
    if (( SHOW_PROGRESS )); then
        if [[ -n $SUBSET_SEQUENCE_LIST ]]; then
            "$IMPG_BIN" similarity -r "$ref:$start-$end" -a "$PAF_FILE" --sequence-files "$FASTA_FILE" --subset-sequence-list "$SUBSET_SEQUENCE_LIST" > "$sim_path"
        else
            "$IMPG_BIN" similarity -r "$ref:$start-$end" -a "$PAF_FILE" --sequence-files "$FASTA_FILE" > "$sim_path"
        fi
    else
        if [[ -n $SUBSET_SEQUENCE_LIST ]]; then
            "$IMPG_BIN" similarity -r "$ref:$start-$end" -a "$PAF_FILE" --sequence-files "$FASTA_FILE" --subset-sequence-list "$SUBSET_SEQUENCE_LIST" > "$sim_path" 2>/dev/null
        else
            "$IMPG_BIN" similarity -r "$ref:$start-$end" -a "$PAF_FILE" --sequence-files "$FASTA_FILE" > "$sim_path" 2>/dev/null
        fi
    fi
done

within_distance() {
    local i=$1
    local j=$2
    local ref1=${REFS[i]}
    local ref2=${REFS[j]}
    if [[ $ref1 != $ref2 ]]; then
        return 1
    fi
    local mid1 mid2 dist
    mid1=$(( (STARTS[i] + ENDS[i]) / 2 ))
    mid2=$(( (STARTS[j] + ENDS[j]) / 2 ))
    if (( mid1 > mid2 )); then
        dist=$(( mid1 - mid2 ))
    else
        dist=$(( mid2 - mid1 ))
    fi
    [[ $dist -le $MAX_DISTANCE ]]
}

run_pair() {
    local i=$1
    local j=$2
    local ref1=${REFS[i]}
    local start1=${STARTS[i]}
    local end1=${ENDS[i]}
    local label1=${LABELS[i]}
    local ref2=${REFS[j]}
    local start2=${STARTS[j]}
    local end2=${ENDS[j]}
    local label2=${LABELS[j]}

    local loc1_sim=${SIM_PATHS[i]}
    local loc2_sim=${SIM_PATHS[j]}
    local pair_dir matrix_path dij_path
    if (( DEBUG_MODE )); then
        pair_dir="$OUT_DIR/${label1}_${label2}"
        mkdir -p "$pair_dir"
        matrix_path="$pair_dir/${label1}_${label2}.haplotype_matrix.tsv"
        dij_path="$pair_dir/${label1}_${label2}.dij.tsv"
    else
        matrix_path=$(mktemp)
        dij_path=$(mktemp)
    fi

    if (( DEBUG_MODE )); then
        log "[$label1/$label2] clustering + haplotype matrix"
        python3 "$CLUSTER_SCRIPT" \
            --locus1 "$loc1_sim" \
            --locus2 "$loc2_sim" \
            --locus1-label "$label1" \
            --locus2-label "$label2" \
            --out-dir "$pair_dir" \
            --haplotype-matrix "$matrix_path" \
            --haplotype-matrix-frequencies
    else
        pair_dir=$(mktemp -d)
        log "[$label1/$label2] clustering + haplotype matrix"
        python3 "$CLUSTER_SCRIPT" \
            --locus1 "$loc1_sim" \
            --locus2 "$loc2_sim" \
            --locus1-label "$label1" \
            --locus2-label "$label2" \
            --out-dir "$pair_dir" \
            --haplotype-matrix "$matrix_path" \
            --haplotype-matrix-frequencies
    fi

    log "[$label1/$label2] computing D_ij + statistics"
    python3 "$DIJ_SCRIPT" --matrix "$matrix_path" --output "$dij_path" --sample-size "$SAMPLE_SIZE" >/dev/null 2>"$TMP_STATS"

    local stat
    stat() {
        local key=$1
        awk -F'\t' -v k="$key" '$1==k {print $2}' "$TMP_STATS"
    }

    local chi_sq cramers_v d_prime r_sq
    chi_sq=$(stat chi_squared)
    cramers_v=$(stat cramers_v)
    d_prime=$(stat composite_d_prime)
    r_sq=$(stat multiallelic_r2)

    if [[ -z $chi_sq || -z $cramers_v || -z $d_prime || -z $r_sq ]]; then
        echo "Error: failed to parse LD statistics for pair ${label1}/${label2}." >&2
        cat "$TMP_STATS" >&2
        exit 1
    fi

    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$ref1" "$start1" "$end1" "$ref2" "$start2" "$end2" "$SAMPLE_SIZE" \
        "$chi_sq" "$cramers_v" "$d_prime" "$r_sq" >> "$SUMMARY_PATH"

    HEATMAP["$i,$j"]=$cramers_v
    HEATMAP["$j,$i"]=$cramers_v

    if (( DEBUG_MODE == 0 )); then
        rm -f "$matrix_path" "$dij_path"
        rm -rf "$pair_dir"
    fi
}

pairs_processed=0
for ((i=0; i<NUM_LOCI; i++)); do
    for ((j=i+1; j<NUM_LOCI; j++)); do
        if within_distance "$i" "$j"; then
            run_pair "$i" "$j"
            pairs_processed=$((pairs_processed + 1))
        else
            log "Skipping ${LABELS[i]} vs ${LABELS[j]}: exceeds $MAX_DISTANCE bp or different contig."
            HEATMAP["$i,$j"]="NA"
            HEATMAP["$j,$i"]="NA"
        fi
    done
done

if [[ $pairs_processed -eq 0 ]]; then
    log "Warning: no locus pairs satisfied the distance/contig criteria."
fi

{
    printf "region"
    for label in "${LABELS[@]}"; do
        printf "\t%s" "$label"
    done
    printf "\n"
    for ((i=0; i<NUM_LOCI; i++)); do
        printf "%s" "${LABELS[i]}"
        for ((j=0; j<NUM_LOCI; j++)); do
            val='.'
            if [[ $j -le $i ]]; then
                val=${HEATMAP["$i,$j"]:-NA}
            fi
            printf "\t%s" "$val"
        done
        printf "\n"
    done
} > "$HEATMAP_PATH"

log "Summary written to $SUMMARY_PATH"
log "Heatmap written to $HEATMAP_PATH"
