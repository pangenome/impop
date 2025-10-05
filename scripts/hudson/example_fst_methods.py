#!/usr/bin/env python3
"""
example_fst_methods.py - Demonstrate the difference between direct and grouped FST calculations
"""

# Create example similarity data
cat > example_similarities.tsv << 'EOF'
group.a	group.b	estimated.identity
seq1_popA	seq2_popA	0.9995
seq1_popA	seq3_popA	0.9993
seq2_popA	seq3_popA	0.9998
seq1_popA	seq4_popB	0.9950
seq1_popA	seq5_popB	0.9948
seq1_popA	seq6_popB	0.9952
seq2_popA	seq4_popB	0.9951
seq2_popA	seq5_popB	0.9949
seq2_popA	seq6_popB	0.9953
seq3_popA	seq4_popB	0.9949
seq3_popA	seq5_popB	0.9947
seq3_popA	seq6_popB	0.9951
seq4_popB	seq5_popB	0.9996
seq4_popB	seq6_popB	0.9994
seq5_popB	seq6_popB	0.9997
EOF

# Create population files
cat > pop_A.txt << 'EOF'
seq1_popA
seq2_popA
seq3_popA
EOF

cat > pop_B.txt << 'EOF'
seq4_popB
seq5_popB
seq6_popB
EOF

echo "Example: Comparing Direct vs Grouped FST Calculations"
echo "===================================================="
echo
echo "Data setup:"
echo "- Population A: 3 sequences (very similar, >0.999 similarity)"
echo "- Population B: 3 sequences (very similar, >0.999 similarity)"
echo "- Between populations: ~0.995 similarity"
echo "- Sequence length: 1,000,000 bp"
echo

echo "1. Direct method (average of all pairwise differences):"
echo "--------------------------------------------------------"
python3 fst.py example_similarities.tsv -a pop_A.txt -b pop_B.txt -l 1000000 -m direct -v
echo

echo "2. Grouped method with threshold 0.999 (groups similar sequences):"
echo "------------------------------------------------------------------"
python3 fst.py example_similarities.tsv -a pop_A.txt -b pop_B.txt -l 1000000 -m grouped -t 0.999 -v
echo

echo "3. Grouped method with threshold 0.996 (less grouping):"
echo "-------------------------------------------------------"
python3 fst.py example_similarities.tsv -a pop_A.txt -b pop_B.txt -l 1000000 -m grouped -t 0.996 -v
echo

echo "Key observations:"
echo "- With threshold 0.999: Sequences within each population are grouped together"
echo "- This effectively treats each population as having a single haplotype"
echo "- The grouped method gives different π values when sequences are collapsed"
echo "- FST values may differ between methods, especially when within-pop diversity is low"
echo
echo "Check the log files for detailed calculations!"

# Running across multiple regions with grouped method
echo
echo "4. Running across multiple regions with grouped method:"
echo "------------------------------------------------------"

cat > test_regions.bed << 'EOF'
chr1	0	1000000
chr2	0	500000
chr3	0	2000000
EOF

echo "Using run_fst.sh with grouped method and threshold 0.999:"
./run_fst.sh -A pop_A.txt -B pop_B.txt -b test_regions.bed -m grouped -t 0.999 -v
