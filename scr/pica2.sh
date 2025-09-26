#!/bin/bash

# Set variables for easier modification
PAF_FILE="../data/hprc465vschm13.aln.paf.gz"
REGION="CHM13#0#chr1:158341439-158343639"
SEQUENCE_FILES="../data/HPRC_r2_assemblies_0.6.1.agc"
SCRIPT_PATH="../scr/pica2.2.py"

# Step 1: Query the GFA and extract window
echo "Step 1: Extracting window from GFA..."
impg query -p "$PAF_FILE" -r "$REGION" --sequence-files "$SEQUENCE_FILES" -o gfa > tmp.1.gfa

# Check if impg query was successful
if [ $? -ne 0 ]; then
    echo "Error: Failed to run impg query"
    exit 1
fi

# Sort and convert the GFA
odgi sort -i tmp.1.gfa -o - | odgi view -i - -g > tmp.gfa

# Check if odgi operations were successful
if [ $? -ne 0 ]; then
    echo "Error: Failed to sort/view GFA"
    exit 1
fi

# Remove intermediate file
rm tmp.1.gfa

# Step 2: Evaluate similarity across pairs of haplotypes
echo "Step 2: Evaluating similarity..."
odgi similarity -i tmp.gfa > tmp.sim

# Check if step 2 was successful
if [ $? -ne 0 ]; then
    echo "Error: Failed to evaluate similarity"
    exit 1
fi

# Step 3: Evaluate nucleotide diversity
echo "Step 3: Calculating nucleotide diversity..."
python3 "$SCRIPT_PATH" tmp.sim -t .988 -l 200 -r 5

# Optional: Clean up temporary files
# rm tmp.gfa tmp.sim

echo "Pipeline completed successfully!"