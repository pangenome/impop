#!/usr/bin/env python3
"""Compute linkage disequilibrium coefficients D_ij and summary metrics."""

import argparse
import csv
import math
import sys


def load_haplotype_matrix(path):
    """Return joint frequencies and marginals parsed from a TSV matrix."""
    with open(path) as fh:
        reader = csv.reader(fh, delimiter='\t')
        header = next(reader)
        if len(header) < 3:
            raise ValueError('Matrix header must have at least two columns plus a total column.')
        col_labels = header[1:-1]

        matrix = {}
        row_totals = {}
        row_labels = []
        col_totals = {label: 0.0 for label in col_labels}

        total_row_encountered = False
        for row in reader:
            label = row[0]
            values = [float(val) for val in row[1:]]
            if label.lower() == 'total':
                if len(values) != len(col_labels) + 1:
                    raise ValueError('Total row must contain one value per column plus the grand total.')
                for col, val in zip(col_labels, values[:-1]):
                    col_totals[col] = val
                grand_total = values[-1]
                total_row_encountered = True
                continue

            if len(values) != len(col_labels) + 1:
                raise ValueError('Each row must contain one value per column plus the marginal total.')

            row_totals[label] = values[-1]
            matrix[label] = {col: val for col, val in zip(col_labels, values[:-1])}
            row_labels.append(label)

        if not total_row_encountered:
            raise ValueError('Input matrix is missing a total row.')

        raw_total = grand_total
        if raw_total == 0:
            raise ValueError('Grand total in matrix cannot be zero.')

        # Normalize to frequencies when totals differ from 1 (counts or rounding drift).
        if abs(raw_total - 1.0) > 1e-9:
            for row in matrix:
                for col in matrix[row]:
                    matrix[row][col] /= raw_total
                row_totals[row] /= raw_total
            for col in col_totals:
                col_totals[col] /= raw_total
            grand_total = 1.0

        return col_labels, row_labels, matrix, row_totals, col_totals, raw_total


def compute_dij(matrix, row_totals, col_totals):
    """Return D_ij values for each row/column pair."""
    result = {}
    for row_label, row_values in matrix.items():
        pi = row_totals[row_label]
        result[row_label] = {}
        for col_label, pij in row_values.items():
            pj = col_totals[col_label]
            result[row_label][col_label] = pij - (pi * pj)
    return result


def write_dij_table(col_labels, row_labels, d_values, path=None):
    """Write D_ij values as a TSV table."""
    output = open(path, 'w', newline='') if path else None
    writer = csv.writer(output or sys.stdout, delimiter='\t')
    writer.writerow(['D_ij'] + col_labels)
    for row_label in row_labels:
        writer.writerow([row_label] + [f"{d_values[row_label][col]:.6f}" for col in col_labels])
    if output:
        output.close()


def compute_chi_square(row_labels, col_labels, d_values, row_totals, col_totals, sample_size):
    component = 0.0
    for row in row_labels:
        pi = row_totals[row]
        if pi == 0:
            continue
        for col in col_labels:
            pj = col_totals[col]
            if pj == 0:
                continue
            dij = d_values[row][col]
            component += (dij ** 2) / (pi * pj)
    return sample_size * component


def compute_cramers_v(chi_sq, sample_size, m, n):
    denom = sample_size * min(m - 1, n - 1)
    if denom <= 0:
        return 0.0
    return math.sqrt(chi_sq / denom)


def d_max(pi, pj, dij):
    if dij >= 0:
        return min(pi * (1 - pj), pj * (1 - pi))
    return min(pi * pj, (1 - pi) * (1 - pj))


def compute_composite_d_prime(row_labels, col_labels, d_values, row_totals, col_totals):
    numerator = 0.0
    denominator = 0.0
    for row in row_labels:
        pi = row_totals[row]
        for col in col_labels:
            pj = col_totals[col]
            dij = d_values[row][col]
            numerator += abs(dij)
            denominator += d_max(pi, pj, dij)
    if denominator == 0:
        return 0.0
    return numerator / denominator


def compute_multiallelic_r2(row_totals, col_totals, d_values, row_labels, col_labels):
    numerator = 0.0
    for row in row_labels:
        for col in col_labels:
            numerator += d_values[row][col] ** 2

    hetero_a = 1 - sum(pi ** 2 for pi in row_totals.values())
    hetero_b = 1 - sum(pj ** 2 for pj in col_totals.values())
    denom = hetero_a * hetero_b
    if denom <= 0:
        return 0.0
    return numerator / denom


def report_statistics(sample_size, chi_sq, cramers_v, d_prime, r_squared):
    print('LD summary statistics', file=sys.stderr)
    print(f"sample_size\t{sample_size:.3f}", file=sys.stderr)
    print(f"chi_squared\t{chi_sq:.6f}", file=sys.stderr)
    print(f"cramers_v\t{cramers_v:.6f}", file=sys.stderr)
    print(f"composite_d_prime\t{d_prime:.6f}", file=sys.stderr)
    print(f"multiallelic_r2\t{r_squared:.6f}", file=sys.stderr)


def parse_args():
    parser = argparse.ArgumentParser(
        description='Compute D_ij = P_ij - p_i * p_j from a haplotype matrix TSV and derive LD summary statistics.'
    )
    parser.add_argument('--matrix', required=True, help='Path to the haplotype frequency matrix TSV.')
    parser.add_argument('--output', help='Optional path to write the D_ij table (stdout if omitted).')
    parser.add_argument(
        '--sample-size',
        type=float,
        help='Sample size N used in chi-square and V computations (defaults to the matrix total).',
    )
    return parser.parse_args()


def main():
    args = parse_args()
    col_labels, row_labels, matrix, row_totals, col_totals, raw_total = load_haplotype_matrix(args.matrix)
    d_values = compute_dij(matrix, row_totals, col_totals)
    write_dij_table(col_labels, row_labels, d_values, args.output)

    if args.sample_size is None and raw_total <= 1.0 + 1e-9:
        print(
            'Warning: sample size defaults to 1. Provide --sample-size when matrix contains normalized frequencies.',
            file=sys.stderr,
        )
    sample_size = args.sample_size if args.sample_size is not None else raw_total
    chi_sq = compute_chi_square(row_labels, col_labels, d_values, row_totals, col_totals, sample_size)
    cramers_v = compute_cramers_v(chi_sq, sample_size, len(row_labels), len(col_labels))
    d_prime = compute_composite_d_prime(row_labels, col_labels, d_values, row_totals, col_totals)
    r_squared = compute_multiallelic_r2(row_totals, col_totals, d_values, row_labels, col_labels)
    report_statistics(sample_size, chi_sq, cramers_v, d_prime, r_squared)


if __name__ == '__main__':
    main()
