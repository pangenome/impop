#!/usr/bin/env python3
import argparse
import csv
import sys
from collections import defaultdict

def load_pairs(path):
    rows = []
    samples = set()
    with open(path) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            a = row['group.a'].split(':', 1)[0]
            b = row['group.b'].split(':', 1)[0]
            val = float(row['estimated.identity'])
            rows.append((a, b, val))
            samples.add(a)
            samples.add(b)
    return rows, sorted(samples)

def union_find(samples):
    parent = {s: s for s in samples}
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    def union(a, b):
        ra, rb = find(a), find(b)
        if ra == rb:
            return
        parent[rb] = ra
    return parent, find, union

def cluster(rows, samples, threshold):
    parent, find, union = union_find(samples)
    for a, b, val in rows:
        if val >= threshold:
            union(a, b)
    comps = defaultdict(list)
    for s in samples:
        comps[find(s)].append(s)
    ordered = sorted(comps.values(), key=lambda c: (-len(c), sorted(c)))
    return ordered

def build_summary(clusters):
    total = sum(len(c) for c in clusters)
    summary = []
    for idx, members in enumerate(clusters, 1):
        cid = f'c{idx}'
        count = len(members)
        freq = (count / total) if total else 0.0
        summary.append((cid, count, freq, sorted(members)))
    return summary

def write_summary(summary, out_file):
    writer = csv.writer(out_file, delimiter='\t')
    writer.writerow(['cluster_id', 'count', 'frequency'])
    for cid, count, freq, _ in summary:
        writer.writerow([cid, count, f"{freq:.6f}"])

def write_details(summary, threshold, path):
    with open(path, 'w', newline='') as fh:
        writer = csv.writer(fh, delimiter='\t')
        writer.writerow(['sample_id', 'cluster_id', 'threshold'])
        for cid, _, _, members in summary:
            for sample in members:
                writer.writerow([sample, cid, threshold])

def main():
    parser = argparse.ArgumentParser(description='Cluster samples in loc.sim-style tables by identity threshold.')
    parser.add_argument('--input', default='loc.sim', help='Path to the similarity table (default: loc.sim)')
    parser.add_argument('--threshold', type=float, default=1.0, help='Minimum estimated.identity to link samples (default: 1.0)')
    parser.add_argument('--output', help='Optional output TSV path for cluster summary; stdout if omitted')
    parser.add_argument('--details', help='Optional path to write detailed sample assignments')
    args = parser.parse_args()

    rows, samples = load_pairs(args.input)
    clusters = cluster(rows, samples, args.threshold)
    summary = build_summary(clusters)

    if args.output:
        with open(args.output, 'w', newline='') as fh:
            write_summary(summary, fh)
    else:
        write_summary(summary, sys.stdout)

    if args.details:
        write_details(summary, args.threshold, args.details)

if __name__ == '__main__':
    main()
