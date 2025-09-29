#!/usr/bin/env python3
"""Plot nucleotide diversity (pi) trends per window for multiple populations."""

import argparse
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import matplotlib.pyplot as plt
import pandas as pd


def chrom_sort_key(chrom: str) -> Tuple[int, str]:
    """Return a sortable key that keeps numeric chromosomes in order."""
    chrom = chrom.replace("chr", "")
    try:
        return (0, int(chrom))
    except ValueError:
        return (1, chrom)


def parse_input(arg: str) -> Tuple[Optional[str], str]:
    """Parse an input argument of the form label=path or just path."""
    if "=" in arg:
        label, path = arg.split("=", 1)
        label = label.strip()
        path = path.strip()
        if not label:
            raise argparse.ArgumentTypeError("Label before '=' must be non-empty")
        return label, path
    return None, arg


def load_pi_table(path: str, label_override: Optional[str]) -> pd.DataFrame:
    """Load a run_pica2_impg output table and return a prepared DataFrame."""
    df = pd.read_csv(path, sep="\t")

    required_columns = {"REGION", "PICA_OUTPUT"}
    missing = required_columns - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns {sorted(missing)} in {path}")

    # Extract chromosome coordinates from REGION values like CHM13#0#chr16:100-200
    region_parts = df["REGION"].str.extract(r"(?:[^#]+#\d+#)?(?P<chrom>[^:]+):(?P<start>\d+)-(?P<end>\d+)")
    if region_parts.isnull().any().any():
        raise ValueError(f"Failed to parse REGION column in {path}")

    df["chrom"] = region_parts["chrom"]
    df["start"] = region_parts["start"].astype(int)
    df["end"] = region_parts["end"].astype(int)
    df["midpoint"] = (df["start"] + df["end"]) / 2.0

    # PICA_OUTPUT field contains e.g. "0.00123 (sequence length: 5000)" – keep the numeric part
    df["pi"] = df["PICA_OUTPUT"].str.split().str[0].astype(float)

    if label_override:
        label = label_override
    elif "SUBSET" in df.columns and df["SUBSET"].nunique() == 1:
        label = str(df["SUBSET"].iloc[0])
    else:
        label = Path(path).stem

    df["label"] = label
    df["source"] = Path(path).as_posix()
    return df


def compute_genome_positions(df: pd.DataFrame) -> pd.DataFrame:
    """Add a genome-wide position for plotting continuous x-axis across chromosomes."""
    ordered_chroms = sorted(df["chrom"].unique(), key=chrom_sort_key)
    offsets: Dict[str, float] = {}
    cumulative = 0.0
    gap = 5e5  # add a gap between chromosomes for visual separation

    for chrom in ordered_chroms:
        chrom_mask = df["chrom"] == chrom
        chrom_end = df.loc[chrom_mask, "end"].max()
        offsets[chrom] = cumulative
        cumulative += chrom_end + gap

    df = df.copy()
    df["genome_pos"] = df["midpoint"] + df["chrom"].map(offsets)
    df["chrom_offset"] = df["chrom"].map(offsets)
    return df


def plot_trend(df: pd.DataFrame, output: Path, title: Optional[str], dpi: int) -> None:
    fig, ax = plt.subplots(figsize=(11, 4))

    for label, subset in df.groupby("label"):
        subset = subset.sort_values("genome_pos")
        ax.plot(
            subset["genome_pos"],
            subset["pi"],
            label=label,
            marker="o",
            markersize=3,
            linewidth=1.2,
        )

    # Chromosome separators and labels
    for chrom, offset in df.drop_duplicates("chrom")[["chrom", "chrom_offset"]].itertuples(index=False):
        ax.axvline(offset, color="grey", linestyle="--", linewidth=0.5, alpha=0.4)

    ax.set_xlabel("Genomic position (window midpoint)")
    ax.set_ylabel("π per window")
    ax.legend(loc="upper right", frameon=False)

    # Use chromosome centers as x ticks
    centers = df.groupby("chrom")["genome_pos"].mean()
    centers = centers.loc[sorted(centers.index, key=chrom_sort_key)]
    ax.set_xticks(centers.values)
    ax.set_xticklabels(centers.index)

    if title:
        ax.set_title(title)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.tight_layout()
    fig.savefig(output, dpi=dpi)
    plt.close(fig)


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Plot π trends from run_pica2_impg output tables.",
    )
    parser.add_argument(
        "-i",
        "--input",
        action="append",
        required=True,
        help="Input file, optionally labelled as label=path. Repeat for multiple populations.",
    )
    parser.add_argument(
        "-o",
        "--output",
        default="pi_trend.png",
        help="Output image path (default: pi_trend.png).",
    )
    parser.add_argument(
        "--title",
        help="Optional plot title.",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=150,
        help="Figure DPI (default: 150).",
    )
    return parser


def main(argv: Optional[List[str]] = None) -> int:
    parser = build_arg_parser()
    args = parser.parse_args(argv)

    records = []
    for item in args.input:
        label, path = parse_input(item)
        df = load_pi_table(path, label)
        records.append(df)

    combined = pd.concat(records, ignore_index=True)
    combined = compute_genome_positions(combined)

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    plot_trend(combined, output_path, args.title, args.dpi)

    print(f"Saved figure to {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
