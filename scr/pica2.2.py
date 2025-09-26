import sys
import argparse
import pandas as pd
import os

def read_similarity_file(filename):
    """
    Read similarity data from file with columns: group.a, group.b, estimated.identity
    Returns list of tuples [(element1, element2, similarity), ...]
    """
    try:
        # Try reading with pandas
        df = pd.read_csv(filename, sep='\t')  # Assuming tab-separated
        
        # Check if required columns exist
        required_cols = ['group.a', 'group.b', 'estimated.identity']
        if not all(col in df.columns for col in required_cols):
            print(f"Error: File must contain columns: {required_cols}")
            print(f"Found columns: {list(df.columns)}")
            sys.exit(1)
        
        # Convert to list of tuples
        similarity_data = []
        for _, row in df.iterrows():
            similarity_data.append((row['group.a'], row['group.b'], row['estimated.identity']))
        
        return similarity_data
        
    except Exception as e:
        print(f"Error reading file {filename}: {e}")
        sys.exit(1)

def analyze_similarity_matrix(similarity_data, threshold=1.0, sequence_length=None, log_file=None, round_digits=None):
    """
    Simple 3-step similarity matrix analysis
    
    Input: 
        similarity_data: list of tuples [(element1, element2, similarity), ...]
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
    
    # Build similarity dictionary
    sim_dict = {}
    elements = set()
    
    for e1, e2, sim in similarity_data:
        # Round similarity if specified
        if round_digits is not None:
            sim = round(sim, round_digits)
            
        sim_dict[(e1, e2)] = sim
        sim_dict[(e2, e1)] = sim  # Make symmetric
        elements.add(e1)
        elements.add(e2)
    
    log_print(f"Loaded {len(similarity_data)} pairwise similarities")
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
            if sim_dict.get((current, other), 0) > threshold:
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
    for i in range(len(groups)):
        for j in range(i+1, len(groups)):
            # Get similarity between groups (any element from each)
            elem1, elem2 = groups[i][0], groups[j][0]
            
            # Check if we have similarity data for this pair
            if (elem1, elem2) in sim_dict:
                similarity = sim_dict[(elem1, elem2)]
            elif (elem2, elem1) in sim_dict:
                similarity = sim_dict[(elem2, elem1)]
            else:
                log_print(f"Warning: No similarity data found between groups G{i+1} and G{j+1}, skipping...")
                continue
            
            # Calculate pair value using (1 - similarity)
            pair_value = (1 - similarity) * len(groups[i]) * len(groups[j])
            group_pairs.append(pair_value)
            log_print(f"  G{i+1}G{j+1}: (1 - {similarity:.6f}) * {len(groups[i])} * {len(groups[j])} = {pair_value:.6f}")
    
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
    similarity_data = read_similarity_file(args.input_file)
    
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
        
        pi, pi_per_site = analyze_similarity_matrix(similarity_data, 
                                                   threshold=args.threshold,
                                                   sequence_length=args.sequence_length,
                                                   log_file=log_file,
                                                   round_digits=args.round_digits)
        
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