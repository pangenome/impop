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


def canonicalize_identifier(identifier: str) -> str:
    """Return a prefix that matches the impg sequence naming scheme.

    The similarity matrix uses identifiers such as:
        HG00097#1#CM094061.1:109468899-109469099

    Population lists contain assembly names such as:
        HG00097_hap1_hprc_r2_v1.0.1
        HG01891_mat_hprc_r2_v1.0.1

    We convert the latter into prefixes (e.g. 'HG00097#1#') that can be used to
    match the similarity records via str.startswith(). When only the sample name
    is available we fall back to 'HG00097#' which matches both haplotypes.
    """

    if not identifier:
        return ""

    token = identifier.strip()
    if not token or token.startswith('#'):
        return ""

    # Remove trailing metadata (e.g. _hprc_r2_v1.0.1)
    if '_hprc' in token:
        token = token.split('_hprc', 1)[0]

    suffix_map = {
        '_hap1': '#1#',
        '_hap2': '#2#',
        '_mat': '#1#',
        '_pat': '#2#',
    }

    for suffix, hap_tag in suffix_map.items():
        if token.endswith(suffix):
            sample = token[:-len(suffix)]
            return f"{sample}{hap_tag}"

    # If the identifier already contains a hap delimiter, keep it as-is
    if '#' in token:
        return token if token.endswith('#') else f"{token}#"

    # Fall back to matching both haplotypes associated with the sample
    return f"{token}#"


def expand_population(raw_ids, all_sequences):
    """Expand population identifiers into the concrete sequence names."""

    expanded = set()
    missing = []

    for raw_id in raw_ids:
        prefix = canonicalize_identifier(raw_id)
        if not prefix:
            continue

        matches = {seq for seq in all_sequences if seq.startswith(prefix)}

        if matches:
            expanded.update(matches)
        else:
            missing.append(raw_id)

    return expanded, missing

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

def calculate_diversity(similarities, seq_set1, seq_set2=None, round_digits=None):
    """
    Calculate average pairwise diversity
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

def calculate_fst(similarities, pop_a, pop_b, sequence_length=None, round_digits=None, log_file=None):
    """Calculate FST using Hudson et al. (1992) formula"""
    
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
    if round_digits is not None:
        log_print(f"Rounding similarities to {round_digits} decimal places")
    log_print("")
    
    # Calculate within-population diversities
    log_print("Within-population diversity (π):")
    pi_a, count_a, miss_a = calculate_diversity(similarities, pop_a, round_digits=round_digits)
    log_print(f"  πA = {pi_a:.6f} (from {count_a} pairs, {miss_a} missing)")
    
    pi_b, count_b, miss_b = calculate_diversity(similarities, pop_b, round_digits=round_digits)
    log_print(f"  πB = {pi_b:.6f} (from {count_b} pairs, {miss_b} missing)")
    
    pi_xy = 0.5 * (pi_a + pi_b)
    log_print(f"  πXY = {pi_xy:.6f} (average of πA and πB)")
    log_print("")
    
    # Calculate between-population diversity
    log_print("Between-population diversity (Dxy):")
    dxy, count_between, miss_between = calculate_diversity(similarities, pop_a, pop_b, round_digits=round_digits)
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
  
Output format:
  FST<tab>pi_A<tab>pi_B<tab>pi_XY<tab>Dxy<tab>Da
  
Where:
  FST = (Dxy - pi_XY) / Dxy  (Hudson et al. 1992)
  pi_A = nucleotide diversity within population A
  pi_B = nucleotide diversity within population B
  pi_XY = (pi_A + pi_B) / 2
  Dxy = nucleotide diversity between populations
  Da = Dxy - pi_XY (net divergence)
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
    pop_a_raw = read_subset_file(args.pop_a)
    pop_b_raw = read_subset_file(args.pop_b)

    pop_a, missing_a = expand_population(pop_a_raw, all_sequences)
    pop_b, missing_b = expand_population(pop_b_raw, all_sequences)

    if args.verbose:
        print(f"Population A candidates: {len(pop_a_raw)}", file=sys.stderr)
        print(f"Population B candidates: {len(pop_b_raw)}", file=sys.stderr)
        print(f"Population A sequences matched: {len(pop_a)}", file=sys.stderr)
        print(f"Population B sequences matched: {len(pop_b)}", file=sys.stderr)

    if missing_a:
        print(
            "Warning: {} identifiers from population A did not match any sequences".format(len(missing_a)),
            file=sys.stderr,
        )
    if missing_b:
        print(
            "Warning: {} identifiers from population B did not match any sequences".format(len(missing_b)),
            file=sys.stderr,
        )

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
            log_file=log_file
        )
    
    # Output results (tab-delimited for easy parsing)
    print(f"{results['fst']:.8f}\t{results['pi_a']:.8f}\t{results['pi_b']:.8f}\t"
          f"{results['pi_xy']:.8f}\t{results['dxy']:.8f}\t{results['da']:.8f}")
    
    if args.verbose:
        print(f"Detailed log saved to: {log_path}", file=sys.stderr)

if __name__ == "__main__":
    main()
