import sys
import argparse
import os
import csv

def read_similarity_file(filename):
    """
    Read similarity data from file with columns: group.a, group.b, estimated.identity
    Returns a tuple (similarity_dict, elements, pair_count) where:
      * similarity_dict maps an unordered pair (min, max) to the similarity value
      * elements is the set of unique group identifiers
      * pair_count counts the number of rows parsed from the file
    """
    try:
        with open(filename, newline='') as handle:
            reader = csv.DictReader(handle, delimiter='\t')

            if reader.fieldnames is None:
                print(f"Error: File {filename} is empty or missing a header")
                sys.exit(1)

            required_cols = {'group.a', 'group.b', 'estimated.identity'}
            missing_cols = required_cols - set(reader.fieldnames)
            if missing_cols:
                print(f"Error: File must contain columns: {sorted(required_cols)}")
                print(f"Found columns: {reader.fieldnames}")
                sys.exit(1)

            similarity_dict = {}
            elements = set()
            pair_count = 0

            for row_number, row in enumerate(reader, start=2):  # start=2 accounts for header row
                pair_count += 1
                e1 = row['group.a']
                e2 = row['group.b']
                try:
                    similarity = float(row['estimated.identity'])
                except (TypeError, ValueError):
                    print(f"Error: Invalid similarity value on line {row_number}: {row['estimated.identity']}")
                    sys.exit(1)

                key = (e1, e2) if e1 <= e2 else (e2, e1)
                similarity_dict[key] = similarity
                elements.add(e1)
                elements.add(e2)

            if pair_count == 0:
                print(f"Warning: No similarity entries found in {filename}")

            return similarity_dict, elements, pair_count

    except FileNotFoundError:
        print(f"Error: File not found {filename}")
        sys.exit(1)
    except Exception as e:
        print(f"Error reading file {filename}: {e}")
        sys.exit(1)

def analyze_similarity_matrix(similarity_dict, elements, pair_count, threshold=1.0, sequence_length=None, log_file=None, round_digits=None):
    """
    Simple 3-step similarity matrix analysis
    
    Input: 
        similarity_dict: mapping of unordered pairs to similarity values
        elements: set of unique element identifiers
        pair_count: number of pairwise similarities parsed from the input file
        threshold: similarity threshold for grouping elements
        sequence_length: length of sequences to normalize pi per site
        log_file: file handle for logging output
        round_digits: number of decimal places to round similarities (None = no rounding)
    Output: pi statistic, pi_per_site
    """
    
    def log_print(message):
        """Print to log file only"""
        if log_file:
            print(message, file=log_file)
    
    # Optionally round similarity values in-place to avoid duplicating data
    if round_digits is not None:
        for key in list(similarity_dict.keys()):
            similarity_dict[key] = round(similarity_dict[key], round_digits)

    def get_similarity(e1, e2):
        key = (e1, e2) if e1 <= e2 else (e2, e1)
        return similarity_dict.get(key)

    log_print(f"Loaded {pair_count} pairwise similarities")
    log_print(f"Found {len(elements)} unique elements")
    if round_digits is not None:
        log_print(f"Rounded similarities to {round_digits} decimal places")
    
    # Step 1: Find groups (elements with similarity > threshold)
    groups = []
    remaining = set(elements)
    
    while remaining:
        # Start new group with first remaining element
        current = remaining.pop()
        group = [current]
        
        # Find all elements similar to current (above threshold)
        for other in list(remaining):
            sim_value = get_similarity(current, other)
            if sim_value is not None and sim_value > threshold:
                group.append(other)
                remaining.remove(other)
        
        groups.append(sorted(group))
    
    groups.sort()  # Sort for consistency
    log_print(f"\nStep 1: Grouping elements (threshold > {threshold})")
    log_print(f"Found {len(groups)} groups:")
    for i, group in enumerate(groups, 1):
        log_print(f"  G{i}: {group} (size: {len(group)})")
    
    # Step 2: Calculate group pairs
    log_print(f"\nStep 2: Calculating group pairs")
    group_pairs = []
    total_elements = sum(len(group) for group in groups)
    if total_elements == 0:
        log_print("Warning: No elements available to compute group pairs")
        return 0.0, 0.0
    for i in range(len(groups)):
        for j in range(i+1, len(groups)):
            # Get similarity between groups (any element from each)
            elem1, elem2 = groups[i][0], groups[j][0]
            
            # Check if we have similarity data for this pair
            similarity = get_similarity(elem1, elem2)
            if similarity is None:
                log_print(f"Warning: No similarity data found between groups G{i+1} and G{j+1}, skipping...")
                continue
            
            # Calculate pair value using (1 - similarity) weighted by group frequencies
            freq_i = len(groups[i]) / total_elements
            freq_j = len(groups[j]) / total_elements
            pair_value = (1 - similarity) * freq_i * freq_j
            group_pairs.append(pair_value)
            log_print(
                f"  G{i+1}G{j+1}: (1 - {similarity:.6f}) * "
                f"({len(groups[i])}/{total_elements}) * "
                f"({len(groups[j])}/{total_elements}) = {pair_value:.6f}"
            )
    
    # Step 3: Calculate pi
    log_print(f"\nStep 3: Calculating pi")
    n = sum(len(group) for group in groups)
    if len(group_pairs) == 0:
        log_print("Warning: No group pairs found with similarity data!")
        return 0.0, 0.0
    
    pi = (n / (n-1)) * sum(2 * pair for pair in group_pairs)
    
    log_print(f"  n (total elements) = {n}")
    log_print(f"  Number of group pairs with data = {len(group_pairs)}")
    log_print(f"  Sum of 2 * group_pairs = {sum(2 * pair for pair in group_pairs):.6f}")
    log_print(f"  pi = {n}/{n-1} * {sum(2 * pair for pair in group_pairs):.6f} = {pi:.6f}")
    
    # Calculate pi per site if sequence length provided
    pi_per_site = None
    if sequence_length:
        pi_per_site = pi / sequence_length
        log_print(f"\nNormalization:")
        log_print(f"  Sequence length = {sequence_length}")
        log_print(f"  pi per site = {pi:.6f} / {sequence_length} = {pi_per_site:.8f}")
    
    return pi, pi_per_site

# Main execution with command line arguments
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Analyze similarity matrix with customizable threshold and sequence length normalization')
    parser.add_argument('input_file', help='Input file with similarity data (TSV format with group.a, group.b, estimated.identity columns)')
    parser.add_argument('--threshold', '-t', type=float, default=0.99, 
                        help='Similarity threshold for grouping elements (default: 0.99)')
    parser.add_argument('--sequence-length', '-l', type=int, 
                        help='Sequence length for normalizing pi per site')
    parser.add_argument('--log-dir', '-d', type=str, default='.',
                        help='Directory to save log file (default: current directory)')
    parser.add_argument('--round-digits', '-r', type=int, default=None,
                        help='Round similarity values to specified decimal places (default: no rounding)')
    
    args = parser.parse_args()
    
    # Create log file name based on input file and specified directory
    base_name = os.path.splitext(os.path.basename(args.input_file))[0]
    log_filename = os.path.join(args.log_dir, f"{base_name}.log")
    
    # Ensure log directory exists
    os.makedirs(args.log_dir, exist_ok=True)
    
    # Read similarity data from file
    similarity_dict, elements, pair_count = read_similarity_file(args.input_file)
    
    # Open log file and run analysis
    with open(log_filename, 'w') as log_file:
        log_file.write(f"Nucleotide Diversity Analysis Log\n")
        log_file.write(f"=================================\n")
        log_file.write(f"Input file: {args.input_file}\n")
        log_file.write(f"Threshold: {args.threshold}\n")
        if args.sequence_length:
            log_file.write(f"Sequence length: {args.sequence_length}\n")
        if args.round_digits is not None:
            log_file.write(f"Similarity rounding: {args.round_digits} decimal places\n")
        log_file.write(f"Log file: {log_filename}\n\n")
        
        pi, pi_per_site = analyze_similarity_matrix(
            similarity_dict,
            elements,
            pair_count,
            threshold=args.threshold,
            sequence_length=args.sequence_length,
            log_file=log_file,
            round_digits=args.round_digits,
        )
        
        log_file.write(f"\n" + "="*50 + "\n")
        log_file.write(f"FINAL RESULTS:\n")
        log_file.write(f"pi = {pi:.6f}\n")
        if pi_per_site is not None:
            log_file.write(f"pi per site = {pi_per_site:.8f}\n")
    
    # Clean console output - only the final result
    if args.sequence_length:
         print(f"{pi_per_site:.8f} (sequence length: {args.sequence_length})")
    else:
        print(f"{pi:.6f} (sequence length: {args.sequence_length})")
    
    # Inform user about log file
   # print(f"# Detailed analysis saved to: {log_filename}", file=sys.stderr)
