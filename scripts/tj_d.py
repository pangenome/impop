#!/usr/bin/env python3
"""
tajimas_d.py — Compute Tajima's D from n (sample size), S (segregating sites), and pi (mean pairwise differences).

Implements the standard formulas (Tajima 1989) using:
    D = (pi - S/a1) / sqrt(e1*S + e2*S*(S-1))
where the constants depend on n:
    a1 = sum_{i=1}^{n-1} 1/i
    a2 = sum_{i=1}^{n-1} 1/i^2
    b1 = (n+1) / [3(n-1)]
    b2 = 2(n^2 + n + 3) / [9n(n-1)]
    c1 = b1 - 1/a1
    c2 = b2 - (n+2)/(a1*n) + a2/(a1**2)
    e1 = c1 / a1
    e2 = c2 / (a1**2 + a2)

Note:
- n must be >= 2.
- S and pi must be >= 0.
- If the denominator is zero (e.g., S == 0), D returns NaN.
"""

from __future__ import annotations
import argparse
import math
from dataclasses import dataclass

@dataclass
class TajimaComponents:
    a1: float
    a2: float
    b1: float
    b2: float
    c1: float
    c2: float
    e1: float
    e2: float
    numerator: float
    denominator: float

def _harmonic(n: int) -> float:
    return sum(1.0 / i for i in range(1, n))

def _harmonic2(n: int) -> float:
    return sum(1.0 / (i * i) for i in range(1, n))

def tajimas_d(n: int, S: float, pi: float, return_components: bool = False):
    if n < 2:
        raise ValueError("n must be >= 2")
    if S < 0 or pi < 0:
        raise ValueError("S and pi must be non-negative")

    a1 = _harmonic(n)
    a2 = _harmonic2(n)
    b1 = (n + 1.0) / (3.0 * (n - 1.0))
    b2 = 2.0 * (n * n + n + 3.0) / (9.0 * n * (n - 1.0))
    c1 = b1 - (1.0 / a1)
    c2 = b2 - ((n + 2.0) / (a1 * n)) + (a2 / (a1 * a1))
    e1 = c1 / a1
    e2 = c2 / (a1 * a1 + a2)

    numerator = pi - (S / a1)
    denominator = math.sqrt(e1 * S + e2 * S * (S - 1.0)) if (S > 0) else float("nan")

    D = numerator / denominator if denominator and not math.isclose(denominator, 0.0) else float("nan")

    if return_components:
        return D, TajimaComponents(a1, a2, b1, b2, c1, c2, e1, e2, numerator, denominator)
    return D

def main():
    parser = argparse.ArgumentParser(description="Compute Tajima's D from n, S, and pi.")
    parser.add_argument("-n", "--sample-size", type=int, required=True, help="Number of sequences (n >= 2)")
    parser.add_argument("-S", "--segregating-sites", type=float, required=True, help="Number of segregating sites S (>= 0)")
    parser.add_argument("-p", "--pi", type=float, required=True, help="Mean pairwise differences pi (>= 0)")
    parser.add_argument("--show-components", action="store_true", help="Print intermediate constants (a1, a2, e1, e2, etc.)")
    args = parser.parse_args()

    D, comps = tajimas_d(args.sample_size, args.segregating_sites, args.pi, return_components=True)

    print(f"Tajima's D: {D}")
    if args.show_components:
        print("--- Components ---")
        print(f"a1={comps.a1} a2={comps.a2}")
        print(f"b1={comps.b1} b2={comps.b2}")
        print(f"c1={comps.c1} c2={comps.c2}")
        print(f"e1={comps.e1} e2={comps.e2}")
        print(f"numerator={comps.numerator} denominator={comps.denominator}")

if __name__ == "__main__":
    main()
