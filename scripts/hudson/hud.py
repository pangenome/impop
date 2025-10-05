#!/usr/bin/env python3
"""
fst.py - Calculate FST (Hudson et al. 1992) from pairwise sequence similarities

FST = (Dxy - πxy) / Dxy

Where:
- Dxy = average pairwise diversity between populations
- πxy = average of within-population diversities
"""

import sys
import argparse
import os
import csv
from collections import defaultdict

def read_similarity_file(filename):
    """Read similarity data from TSV file"""
    try:
        with open(filename, newline='') as f:
            reader = csv.DictReader(f, delimiter='\t')
            
            if not reader.fieldnames:
                print(f"Error: Empty file {filename}", file=sys.stderr)
                sys.exit(1)
                
            required = {'group.a', 'group.b', 'estimated.identity'}
            if not required.issubset(set(reader.fieldnames)):
                print(f"Error: File must contain columns: {required}", file=sys.stderr)
                print(f"Found: {reader.fieldnames}", file=sys.stderr)
                sys.exit(1)
            
            similarities = {}
            all_sequences = set()
            
            for row in reader:
                seq1, seq2 = row['group.a'], row['group.b']
                try:
                    sim = float(row['estimated.identity'])
                except ValueError:
                    print(f"Warning: Invalid similarity value: {row['estimated.identity']}", file=sys.stderr)
                    continue
                
                key = (seq1, seq2) if seq1 <= seq2 else (seq2, seq1)
                similarities[key] = sim
                all_sequences.update([seq1, seq2])
            
            return similarities, all_sequences
            
    except FileNotFoundError:
        print(f"Error: File not found: {filename}", file=sys.stderr)
        sys.exit(1)

def read_subset_file(filename):
    """Read sequence IDs from a file"""
    try:
        with open(filename) as f:
            return set(line.strip() for line in f if line.strip() and not line.startswith('#'))
    except FileNotFoundError:
        print(f"Error: Subset file not found: {filename}", file=sys.stderr)
        sys.exit(1)

def group_sequences(similarities, sequences, threshold=0.999, round_digits=None):
    """Group sequences with similarity > threshold"""
    groups = []
    remaining = set(sequences)
    
    while remaining:
        current = remaining.pop()
        group = [current]
        
        # Find all sequences similar to current
        for other in list(remaining):
            key = (current, other) if current <= other else (other, current)
            if key in similarities:
                sim = similarities[key]
                if round_digits is not None:
                    sim = round(sim, round_digits)
                if sim > threshold:
                    group.append(other)
                    remaining.remove(other)
        
        groups.append(sorted(group))
    
    return sorted(groups)

def get_group_similarity(similarities, group1, group2, round_digits=None):
    """Get representative similarity between two groups (using first found pair)"""
    for seq1 in group1:
        for seq2 in group2:
            key = (seq1, seq2) if seq1 <= seq2 else (seq2, seq1)
            if key in similarities:
                sim = similarities[key]
                if round_digits is not None:
                    sim = round(sim, round_digits)
                return sim
    return None

def calculate_diversity_grouped(similarities, sequences, threshold=0.999, round_digits=None):
    """Calculate diversity using frequency-based formula after grouping"""
    groups = group_sequences(similarities, sequences, threshold, round_digits)
    n_total = len(sequences)
    
    if n_total <= 1:
        return 0.0, len(groups), 0
    
    # Calculate diversity between groups
    diversity_sum = 0.0
    pair_count = 0
    missing_count = 0
    
    for i in range(len(groups)):
        for j in range(i + 1, len(groups)):
            sim = get_group_similarity(similarities, groups[i], groups[j], round_digits)
            
            if sim is not None:
                freq_i = len(groups[i]) / n_total
                freq_j = len(groups[j]) / n_total
                diversity_sum += 2 * freq_i * freq_j * (1 - sim)
                pair_count += 1
            else:
                missing_count += 1
    
    # Apply Bessel correction
    diversity = diversity_sum * n_total / (n_total - 1)
    
    return diversity, len(groups), missing_count

def calculate_diversity_direct(similarities, seq_set1, seq_set2=None, round_digits=None):
    """
    Calculate average pairwise diversity (direct method)
    If seq_set2 is None: calculate within seq_set1
    If seq_set2 is provided: calculate between seq_set1 and seq_set2
    """
    diversities = []
    missing = 0
    
    if seq_set2 is None:
        # Within-population diversity
        seq_list = list(seq_set1)
        for i in range(len(seq_list)):
            for j in range(i + 1, len(seq_list)):
                seq1, seq2 = seq_list[i], seq_list[j]
                key = (seq1, seq2) if seq1 <= seq2 else (seq2, seq1)
                
                if key in similarities:
                    sim = similarities[key]
                    if round_digits is not None:
                        sim = round(sim, round_digits)
                    diversities.append(1 - sim)
                else:
                    missing += 1
    else:
        # Between-population diversity
        for seq1 in seq_set1:
            for seq2 in seq_set2:
                key = (seq1, seq2) if seq1 <= seq2 else (seq2, seq1)
                
                if key in similarities:
                    sim = similarities[key]
                    if round_digits is not None:
                        sim = round(sim, round_digits)
                    diversities.append(1 - sim)
                else:
                    missing += 1
    
    if not diversities:
        return 0.0, 0, missing
    
    return sum(diversities) / len(diversities), len(diversities), missing

def calculate_fst(similarities, pop_a, pop_b, sequence_length=None, round_digits=None, 
                  log_file=None, method='direct', threshold=0.999):
    """
    Calculate FST using Hudson et al. (1992) formula
    
    method: 'direct' for pairwise average, 'grouped' for frequency-based with grouping
    threshold: similarity threshold for grouping (only used when method='grouped')
    """
    
    def log_print(msg):
        if log_file:
            print(msg, file=log_file)
    
    # Ensure populations don't overlap
    overlap = pop_a & pop_b
    if overlap:
        print(f"Warning: {len(overlap)} sequences appear in both populations", file=sys.stderr)
        pop_a = pop_a - overlap
        pop_b = pop_b - overlap
    
    log_print("FST Calculation")
    log_print("=" * 50)
    log_print(f"Population A: {len(pop_a)} sequences")
    log_print(f"Population B: {len(pop_b)} sequences")
    log_print(f"Method: {method}")
    if method == 'grouped':
        log_print(f"Grouping threshold: {threshold}")
    if round_digits is not None:
        log_print(f"Rounding similarities to {round_digits} decimal places")
    log_print("")
    
    # Calculate within-population diversities
    if method == 'grouped':
        log_print("Within-population diversity (π) using grouped method:")
        pi_a, groups_a, miss_a = calculate_diversity_grouped(
            similarities, pop_a, threshold, round_digits
        )
        log_print(f"  πA = {pi_a:.6f} ({groups_a} groups from {len(pop_a)} sequences, {miss_a} missing pairs)")
        
        pi_b, groups_b, miss_b = calculate_diversity_grouped(
            similarities, pop_b, threshold, round_digits
        )
        log_print(f"  πB = {pi_b:.6f} ({groups_b} groups from {len(pop_b)} sequences, {miss_b} missing pairs)")
    else:
        log_print("Within-population diversity (π) using direct method:")
        pi_a, count_a, miss_a = calculate_diversity_direct(
            similarities, pop_a, round_digits=round_digits
        )
        log_print(f"  πA = {pi_a:.6f} (from {count_a} pairs, {miss_a} missing)")
        
        pi_b, count_b, miss_b = calculate_diversity_direct(
            similarities, pop_b, round_digits=round_digits
        )
        log_print(f"  πB = {pi_b:.6f} (from {count_b} pairs, {miss_b} missing)")
    
    pi_xy = 0.5 * (pi_a + pi_b)
    log_print(f"  πXY = {pi_xy:.6f} (average of πA and πB)")
    log_print("")
    
    # Calculate between-population diversity
    log_print("Between-population diversity (Dxy):")
    
    if method == 'grouped':
        # For Dxy with grouping, we need to handle between-population groups specially
        groups_a = group_sequences(similarities, pop_a, threshold, round_digits)
        groups_b = group_sequences(similarities, pop_b, threshold, round_digits)
        
        n_total_a = len(pop_a)
        n_total_b = len(pop_b)
        n_total = n_total_a + n_total_b
        
        dxy_sum = 0.0
        pair_count = 0
        missing_count = 0
        
        for group_a in groups_a:
            for group_b in groups_b:
                sim = get_group_similarity(similarities, group_a, group_b, round_digits)
                if sim is not None:
                    freq_a = len(group_a) / n_total
                    freq_b = len(group_b) / n_total
                    # For between-population diversity, we use the product of frequencies
                    # normalized by the relative population sizes
                    weight = (len(group_a) * len(group_b)) / (n_total_a * n_total_b)
                    dxy_sum += weight * (1 - sim)
                    pair_count += 1
                else:
                    missing_count += 1
        
        dxy = dxy_sum
        log_print(f"  Dxy = {dxy:.6f} (from {len(groups_a)} x {len(groups_b)} group pairs, {missing_count} missing)")
    else:
        dxy, count_between, miss_between = calculate_diversity_direct(
            similarities, pop_a, pop_b, round_digits=round_digits
        )
        log_print(f"  Dxy = {dxy:.6f} (from {count_between} pairs, {miss_between} missing)")
    
    log_print("")
    
    # Calculate FST
    if dxy > 0:
        fst = (dxy - pi_xy) / dxy
        log_print("FST calculation:")
        log_print(f"  FST = (Dxy - πXY) / Dxy")
        log_print(f"      = ({dxy:.6f} - {pi_xy:.6f}) / {dxy:.6f}")
        log_print(f"      = {fst:.6f}")
    else:
        fst = 0.0
        log_print("FST = 0 (Dxy = 0)")
    
    # Per-site calculations if sequence length provided
    if sequence_length and sequence_length > 0:
        log_print("")
        log_print(f"Per-site values (sequence length = {sequence_length:,}):")
        log_print(f"  πA per site = {pi_a/sequence_length:.8f}")
        log_print(f"  πB per site = {pi_b/sequence_length:.8f}")
        log_print(f"  πXY per site = {pi_xy/sequence_length:.8f}")
        log_print(f"  Dxy per site = {dxy/sequence_length:.8f}")
        
        return {
            'fst': fst,
            'pi_a': pi_a / sequence_length,
            'pi_b': pi_b / sequence_length,
            'pi_xy': pi_xy / sequence_length,
            'dxy': dxy / sequence_length,
            'da': (dxy - pi_xy) / sequence_length
        }
    else:
        return {
            'fst': fst,
            'pi_a': pi_a,
            'pi_b': pi_b,
            'pi_xy': pi_xy,
            'dxy': dxy,
            'da': dxy - pi_xy
        }

def main():
    parser = argparse.ArgumentParser(
        description='Calculate FST from pairwise sequence similarities',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
  %(prog)s similarities.tsv -a pop_a.txt -b pop_b.txt -l 1000000
  %(prog)s similarities.tsv -a pop_a.txt -b pop_b.txt -l 1000000 -m grouped -t 0.999
  
Output format:
  FST<tab>pi_A<tab>pi_B<tab>pi_XY<tab>Dxy<tab>Da
  
Where:
  FST = (Dxy - pi_XY) / Dxy  (Hudson et al. 1992)
  pi_A = nucleotide diversity within population A
  pi_B = nucleotide diversity within population B
  pi_XY = (pi_A + pi_B) / 2
  Dxy = nucleotide diversity between populations
  Da = Dxy - pi_XY (net divergence)
  
Methods:
  direct: Average pairwise diversity (default)
  grouped: Frequency-based calculation after grouping similar sequences
        """
    )
    
    parser.add_argument('similarity_file', 
                        help='TSV file with columns: group.a, group.b, estimated.identity')
    parser.add_argument('-a', '--pop-a', required=True,
                        help='File listing sequence IDs for population A')
    parser.add_argument('-b', '--pop-b', required=True,
                        help='File listing sequence IDs for population B')
    parser.add_argument('-l', '--length', type=int, default=None,
                        help='Sequence length for per-site calculations')
    parser.add_argument('-r', '--round', type=int, default=None,
                        help='Round similarities to N decimal places')
    parser.add_argument('-m', '--method', choices=['direct', 'grouped'], default='direct',
                        help='Calculation method: direct or grouped (default: direct)')
    parser.add_argument('-t', '--threshold', type=float, default=0.999,
                        help='Similarity threshold for grouping (default: 0.999, used only with -m grouped)')
    parser.add_argument('-d', '--log-dir', default='.',
                        help='Directory for log file (default: current directory)')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='Print detailed progress to stderr')
    
    args = parser.parse_args()
    
    # Read input files
    if args.verbose:
        print(f"Reading similarity file: {args.similarity_file}", file=sys.stderr)
    similarities, all_sequences = read_similarity_file(args.similarity_file)
    
    if args.verbose:
        print(f"Reading population files...", file=sys.stderr)
    pop_a = read_subset_file(args.pop_a)
    pop_b = read_subset_file(args.pop_b)
    
    if args.verbose:
        print(f"Population A: {len(pop_a)} sequences", file=sys.stderr)
        print(f"Population B: {len(pop_b)} sequences", file=sys.stderr)
        print(f"Method: {args.method}", file=sys.stderr)
        if args.method == 'grouped':
            print(f"Grouping threshold: {args.threshold}", file=sys.stderr)
    
    # Verify sequences exist in similarity file
    missing_a = pop_a - all_sequences
    missing_b = pop_b - all_sequences
    
    if missing_a:
        print(f"Warning: {len(missing_a)} sequences from population A not found in similarity file", 
              file=sys.stderr)
    if missing_b:
        print(f"Warning: {len(missing_b)} sequences from population B not found in similarity file", 
              file=sys.stderr)
    
    # Remove missing sequences
    pop_a = pop_a & all_sequences
    pop_b = pop_b & all_sequences
    
    if not pop_a or not pop_b:
        print("Error: No valid sequences found in one or both populations", file=sys.stderr)
        sys.exit(1)
    
    # Create log file
    base_name = os.path.splitext(os.path.basename(args.similarity_file))[0]
    log_path = os.path.join(args.log_dir, f"{base_name}_fst.log")
    os.makedirs(args.log_dir, exist_ok=True)
    
    # Calculate FST
    with open(log_path, 'w') as log_file:
        results = calculate_fst(
            similarities, pop_a, pop_b, 
            sequence_length=args.length,
            round_digits=args.round,
            log_file=log_file,
            method=args.method,
            threshold=args.threshold
        )
    
    # Output results (tab-delimited for easy parsing)
    print(f"{results['fst']:.8f}\t{results['pi_a']:.8f}\t{results['pi_b']:.8f}\t"
          f"{results['pi_xy']:.8f}\t{results['dxy']:.8f}\t{results['da']:.8f}")
    
    if args.verbose:
        print(f"Detailed log saved to: {log_path}", file=sys.stderr)

if __name__ == "__main__":
    main()