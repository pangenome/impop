#!/bin/bash

# Set variables for easier modification
BED_FILE="$1"  # Take BED file as first argument
PAF_FILE="../data/hprc465vschm13.aln.paf.gz"
SEQUENCE_FILES="../data/HPRC_r2_assemblies_0.6.1.agc"
SCRIPT_PATH="../impop/scr/pica2.2.py"

# Check if BED file is provided
if [ -z "$1" ]; then
    echo "Usage: $0 <bed_file>"
    echo "Example: $0 ackr1_win.bed"
    exit 1
fi

# Check if BED file exists
if [ ! -f "$BED_FILE" ]; then
    echo "Error: BED file '$BED_FILE' not found"
    exit 1
fi

# Print header
echo -e "REGION\tPICA_OUTPUT"

# Process each line in the BED file
while IFS=$'\t' read -r chr start end name; do
    # Skip header lines or empty lines
    if [[ "$chr" == "#"* ]] || [ -z "$chr" ]; then
        continue
    fi
    
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
    
    # Step 3: Evaluate nucleotide diversity and capture output
    PICA_OUTPUT=$(python3 "$SCRIPT_PATH" tmp.sim -t .988 -l 200 -r 5 2>&1 | tr '\n' ' ' | sed 's/  */ /g' | sed 's/ *$//')
    
    # Print region and pica output on same line
    echo -e "${REGION}\t${PICA_OUTPUT}"
    
    # Clean up temporary files before processing next region
    rm -f tmp.gfa tmp.sim
    
done < "$BED_FILE"