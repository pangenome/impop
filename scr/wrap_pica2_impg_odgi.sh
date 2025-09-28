#!/bin/bash

# Default values
PAF_FILE="../data/hprc465vschm13.aln.paf.gz"
SEQUENCE_FILES="../data/HPRC_r2_assemblies_0.6.1.agc"
SCRIPT_PATH="pica2.2.py"

# Function to display usage
usage() {
    echo "Usage: $0 -b <bed_file> -t <threshold> -r <r_value>"
    echo "  -b  BED file containing regions (required)"
    echo "  -t  Similarity threshold for pica2.py (required)"
    echo "  -r  R value for pica2.py (required)"
    echo "  -h  Display this help message"
    echo ""
    echo "Example: $0 -b ackr1_win.bed -t 0.988 -r 5"
    exit 1
}

# Parse command line arguments
while getopts "b:t:r:h" opt; do
    case $opt in
        b) BED_FILE="$OPTARG" ;;
        t) THRESHOLD="$OPTARG" ;;
        r) R_VALUE="$OPTARG" ;;
        h) usage ;;
        *) usage ;;
    esac
done

# Check if all required arguments are provided
if [ -z "$BED_FILE" ] || [ -z "$THRESHOLD" ] || [ -z "$R_VALUE" ]; then
    echo "Error: Missing required arguments"
    usage
fi

# Check if BED file exists
if [ ! -f "$BED_FILE" ]; then
    echo "Error: BED file '$BED_FILE' not found"
    exit 1
fi

# Validate threshold is a number
if ! [[ "$THRESHOLD" =~ ^[0-9]*\.?[0-9]+$ ]]; then
    echo "Error: Threshold must be a number (e.g., 0.988)"
    exit 1
fi

# Validate r value is a number
if ! [[ "$R_VALUE" =~ ^[0-9]+$ ]]; then
    echo "Error: R value must be an integer (e.g., 5)"
    exit 1
fi

# Print header
echo -e "REGION\tLENGTH\tTHRESHOLD\tR_VALUE\tPICA_OUTPUT"

# Process each line in the BED file
while IFS=$'\t' read -r chr start end name; do
    # Skip header lines or empty lines
    if [[ "$chr" == "#"* ]] || [ -z "$chr" ]; then
        continue
    fi
    
    # Calculate window length
    LENGTH=$((end - start))
    
    # Construct the region string
    REGION="CHM13#0#${chr}:${start}-${end}"
    
    # Suppress verbose output by redirecting to stderr
    {
        # Step 1: Query the GFA and extract window
        impg query -p "$PAF_FILE" -r "$REGION" --sequence-files "$SEQUENCE_FILES" -o gfa > tmp.1.gfa
        
        # Check if impg query was successful
        if [ $? -ne 0 ]; then
            echo "Error: Failed to run impg query for region $REGION" >&2
            continue
        fi
        
        # Sort and convert the GFA
        odgi sort -i tmp.1.gfa -o - | odgi view -i - -g > tmp.gfa
        
        # Check if odgi operations were successful
        if [ $? -ne 0 ]; then
            echo "Error: Failed to sort/view GFA for region $REGION" >&2
            rm -f tmp.1.gfa
            continue
        fi
        
        # Remove intermediate file
        rm tmp.1.gfa
        
        # Step 2: Evaluate similarity across pairs of haplotypes
        odgi similarity -i tmp.gfa > tmp.sim
        
        # Check if step 2 was successful
        if [ $? -ne 0 ]; then
            echo "Error: Failed to evaluate similarity for region $REGION" >&2
            rm -f tmp.gfa
            continue
        fi
    } 2>/dev/null  # Suppress standard error for cleaner output
    
    # Step 3: Evaluate nucleotide diversity with calculated length
    PICA_OUTPUT=$(python3 "$SCRIPT_PATH" tmp.sim -t "$THRESHOLD" -l "$LENGTH" -r "$R_VALUE" 2>&1 | tr '\n' ' ' | sed 's/  */ /g' | sed 's/ *$//')
    
    # Print region, parameters, and pica output on same line
    echo -e "${REGION}\t${LENGTH}\t${THRESHOLD}\t${R_VALUE}\t${PICA_OUTPUT}"
    
    # Clean up temporary files before processing next region
    rm -f tmp.gfa tmp.sim
    
done < "$BED_FILE"