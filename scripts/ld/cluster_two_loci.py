#!/usr/bin/env python3
"""Compute per-locus clusters, merged genotypes, and haplotype frequencies."""

import argparse
import csv
import os
from collections import Counter

import cluster_samples as cs


def apply_prefix(summary, prefix):
    """Return a summary with allele IDs prefixed (e.g. l1c1)."""
    prefixed = []
    for cid, count, freq, members in summary:
        new_cid = f"{prefix}{cid}"
        prefixed.append((new_cid, count, freq, members))
    return prefixed


def cluster_locus(path, threshold, prefix):
    """Return cluster summary and sample assignments for a single locus."""
    rows, samples = cs.load_pairs(path)
    clusters = cs.cluster(rows, samples, threshold)
    summary = cs.build_summary(clusters)
    summary = apply_prefix(summary, prefix)
    assignments = {}
    for cid, _, _, members in summary:
        for sample in members:
            assignments[sample] = cid
    return summary, assignments


def write_summary_table(summary, path):
    """Write the allele frequency table for a locus."""
    with open(path, 'w', newline='') as fh:
        cs.write_summary(summary, fh)


def write_genotype_table(summary, threshold, path):
    """Write the per-sample genotype assignments for a locus."""
    cs.write_details(summary, threshold, path)


def write_merged_genotypes(assignments1, assignments2, label1, label2, path):
    """Write merged genotypes across two loci for all observed samples."""
    samples = sorted(set(assignments1) & set(assignments2))
    with open(path, 'w', newline='') as fh:
        writer = csv.writer(fh, delimiter='\t')
        writer.writerow(['sample_id', label1, label2])
        for sample in samples:
            call1 = assignments1.get(sample, 'NA')
            call2 = assignments2.get(sample, 'NA')
            writer.writerow([sample, call1, call2])


def collect_haplotypes(assignments1, assignments2):
    """Return haplotype counts for samples present at both loci."""
    samples = set(assignments1) & set(assignments2)
    counter = Counter()
    for sample in samples:
        counter[(assignments1[sample], assignments2[sample])] += 1
    return counter


def write_haplotype_frequencies(counter, label1, label2, path):
    """Write haplotype frequencies computed from merged genotype calls."""
    with open(path, 'w', newline='') as fh:
        writer = csv.writer(fh, delimiter='\t')
        writer.writerow([f'{label1}{label2}', 'count'])
        for (allele1, allele2), count in sorted(counter.items(), key=lambda item: (-item[1], item[0][0], item[0][1])):
            writer.writerow([f"{allele1}{allele2}", count])


def write_haplotype_matrix(counter, summary1, summary2, label1, label2, path, use_frequencies=False):
    """Write a contingency matrix of allele co-occurrence with marginal sums."""
    row_labels = [cid for cid, *_ in summary1]
    col_labels = [cid for cid, *_ in summary2]

    col_totals = {col: 0 for col in col_labels}
    grand_total = sum(counter.values())

    def format_value(value):
        if not use_frequencies:
            return value
        if grand_total == 0:
            return "0.000000"
        return f"{(value / grand_total):.6f}"

    with open(path, 'w', newline='') as fh:
        writer = csv.writer(fh, delimiter='\t')
        writer.writerow([label1] + col_labels + ['total'])

        for row in row_labels:
            row_counts = []
            row_total = 0
            for col in col_labels:
                value = counter.get((row, col), 0)
                row_counts.append(value)
                row_total += value
                col_totals[col] += value
            writer.writerow([row] + [format_value(v) for v in row_counts] + [format_value(row_total)])

        total_row = [format_value(col_totals[col]) for col in col_labels]
        writer.writerow(['total'] + total_row + [format_value(grand_total)])


def parse_args():
    parser = argparse.ArgumentParser(
        description='Compute allele, genotype, merged genotype, and haplotype summaries for two loci.'
    )
    parser.add_argument('--locus1', required=True, help='Path to the first locus similarity table (TSV).')
    parser.add_argument('--locus2', required=True, help='Path to the second locus similarity table (TSV).')
    parser.add_argument('--locus1-label', default='loc1', help='Label used for the first locus in outputs.')
    parser.add_argument('--locus2-label', default='loc2', help='Label used for the second locus in outputs.')
    parser.add_argument('--locus1-allele-prefix', default='l1', help='Prefix applied to allele IDs for the first locus (default: l1).')
    parser.add_argument('--locus2-allele-prefix', default='l2', help='Prefix applied to allele IDs for the second locus (default: l2).')
    parser.add_argument('--threshold', type=float, default=1.0, help='Default identity threshold for clustering (default: 1.0).')
    parser.add_argument('--locus1-threshold', type=float, help='Optional locus-specific threshold (overrides --threshold).')
    parser.add_argument('--locus2-threshold', type=float, help='Optional locus-specific threshold (overrides --threshold).')
    parser.add_argument('--out-dir', default='.', help='Directory where output tables are written.')
    parser.add_argument('--merged-genotypes', help='Optional path for merged genotypes table.')
    parser.add_argument('--haplotype-frequencies', help='Optional path for haplotype frequency table.')
    parser.add_argument('--haplotype-matrix', help='Optional path for haplotype contingency matrix.')
    parser.add_argument(
        '--haplotype-matrix-frequencies',
        action='store_true',
        help='If set, write haplotype matrix values as frequencies instead of counts.',
    )
    return parser.parse_args()


def main():
    args = parse_args()

    locus1_threshold = args.locus1_threshold if args.locus1_threshold is not None else args.threshold
    locus2_threshold = args.locus2_threshold if args.locus2_threshold is not None else args.threshold

    summary1, assignments1 = cluster_locus(args.locus1, locus1_threshold, args.locus1_allele_prefix)
    summary2, assignments2 = cluster_locus(args.locus2, locus2_threshold, args.locus2_allele_prefix)

    os.makedirs(args.out_dir, exist_ok=True)

    freq1_path = os.path.join(args.out_dir, f"{args.locus1_label}.freq")
    freq2_path = os.path.join(args.out_dir, f"{args.locus2_label}.freq")
    geno1_path = os.path.join(args.out_dir, f"{args.locus1_label}.geno")
    geno2_path = os.path.join(args.out_dir, f"{args.locus2_label}.geno")
    merged_path = args.merged_genotypes or os.path.join(args.out_dir, f"{args.locus1_label}_{args.locus2_label}.merged.tsv")
    haplo_path = args.haplotype_frequencies or os.path.join(args.out_dir, f"{args.locus1_label}{args.locus2_label}.haplotypes.tsv")

    write_summary_table(summary1, freq1_path)
    write_summary_table(summary2, freq2_path)
    write_genotype_table(summary1, locus1_threshold, geno1_path)
    write_genotype_table(summary2, locus2_threshold, geno2_path)
    write_merged_genotypes(assignments1, assignments2, args.locus1_label, args.locus2_label, merged_path)

    hap_counter = collect_haplotypes(assignments1, assignments2)
    write_haplotype_frequencies(hap_counter, args.locus1_label, args.locus2_label, haplo_path)

    matrix_path = args.haplotype_matrix or os.path.join(
        args.out_dir,
        f"{args.locus1_label}_{args.locus2_label}.haplotype_matrix.tsv",
    )
    write_haplotype_matrix(
        hap_counter,
        summary1,
        summary2,
        args.locus1_label,
        args.locus2_label,
        matrix_path,
        use_frequencies=args.haplotype_matrix_frequencies,
    )


if __name__ == '__main__':
    main()
