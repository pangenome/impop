# How linkage disequilibrium is evaluated

This document describes the current LD workflow implemented in
[`scripts/ld/run_ld_pipeline.sh`](/home/enza/iusatoProtocol/imscratch/impop/scripts/ld/run_ld_pipeline.sh),
[`scripts/ld/cluster_two_loci.py`](/home/enza/iusatoProtocol/imscratch/impop/scripts/ld/cluster_two_loci.py),
and [`scripts/ld/compute_dij.py`](/home/enza/iusatoProtocol/imscratch/impop/scripts/ld/compute_dij.py).

## Mathematical definition

For two multiallelic loci:

- locus A has alleles \(A_1, A_2, \ldots, A_m\)
- locus B has alleles \(B_1, B_2, \ldots, B_n\)

Let \(P_{ij}\) be the joint haplotype frequency of allele pair \((A_i, B_j)\).

The joint frequencies satisfy:

\[
\sum_{i=1}^{m}\sum_{j=1}^{n} P_{ij} = 1
\]

Marginal allele frequencies are:

\[
p_{i\cdot} = \sum_{j=1}^{n} P_{ij}
\]

\[
p_{\cdot j} = \sum_{i=1}^{m} P_{ij}
\]

The disequilibrium coefficient for each allele pair is:

\[
D_{ij} = P_{ij} - p_{i\cdot} p_{\cdot j}
\]

## Summary statistics

Because multiallelic loci produce many \(D_{ij}\) values, the workflow reports summary measures across the whole contingency table.

### Chi-square

\[
\chi^2 = N \sum_{i=1}^{m}\sum_{j=1}^{n} \frac{D_{ij}^2}{p_{i\cdot} p_{\cdot j}}
\]

where \(N\) is the sample size used for the LD summary.

### Cramer's V

\[
V = \sqrt{\frac{\chi^2}{N \cdot \min(m - 1, n - 1)}}
\]

### Composite \(D'\)

The implementation reports a multiallelic composite \(D'\):

\[
D' = \frac{\sum_{i=1}^{m}\sum_{j=1}^{n} |D_{ij}|}
{\sum_{i=1}^{m}\sum_{j=1}^{n} D^{\max}_{ij}}
\]

where \(D^{\max}_{ij}\) depends on the sign of \(D_{ij}\):

- if \(D_{ij} \ge 0\), \(D^{\max}_{ij} = \min(p_{i\cdot}(1-p_{\cdot j}), p_{\cdot j}(1-p_{i\cdot}))\)
- if \(D_{ij} < 0\), \(D^{\max}_{ij} = \min(p_{i\cdot}p_{\cdot j}, (1-p_{i\cdot})(1-p_{\cdot j}))\)

### Multiallelic \(r^2\)

The implementation uses generalized heterozygosity terms:

\[
H_A = 1 - \sum_i p_{i\cdot}^2
\]

\[
H_B = 1 - \sum_j p_{\cdot j}^2
\]

\[
r^2 = \frac{\sum_{i=1}^{m}\sum_{j=1}^{n} D_{ij}^2}{H_A H_B}
\]

## Current workflow

The recommended entrypoint is:

[`scripts/ld/run_ld_pipeline.sh`](/home/enza/iusatoProtocol/imscratch/impop/scripts/ld/run_ld_pipeline.sh)

The pipeline does this:

1. Generate one `impg similarity` table per locus.
2. For each eligible locus pair, cluster alleles and build a normalized haplotype matrix.
3. Compute `D_ij`, chi-square, Cramer's V, composite `D'`, and multiallelic `r^2`.
4. Write per-pair summary rows to `ld_summary.tsv` and the lower-triangular Cramer's V matrix to `ld_heatmap.tsv`.

By default the script is silent and keeps only:

- `ld_summary.tsv`
- `ld_heatmap.tsv`
- `loci/*.sim`

If `--debug` is used, it also keeps per-pair intermediate directories such as `loc1_loc2/`.

## Input modes

The wrapper supports two mutually exclusive region input modes.

### 1. BED mode

```bash
bash scripts/ld/run_ld_pipeline.sh \
  --bed regions.bed \
  -p ../data/hprc465vschm13.aln.paf.gz \
  -f ../data/HPRC_r2_assemblies_0.6.1.agc \
  -n 233
```

### 2. Tiled region mode

```bash
bash scripts/ld/run_ld_pipeline.sh \
  --region CHM13#0#chr22:11582486-11589486 \
  --interval-size 1000 \
  -p ../data/hprc465vschm13.aln.paf.gz \
  -f ../data/HPRC_r2_assemblies_0.6.1.agc \
  -n 233
```

This splits the requested range into consecutive intervals of `--interval-size` bp and labels them `loc1`, `loc2`, and so on.

## Optional subsetting

If only a subset of sequences should be considered, pass a file containing sequence names:

```bash
bash scripts/ld/run_ld_pipeline.sh \
  --bed regions.bed \
  -p ../data/hprc465vschm13.aln.paf.gz \
  -f ../data/HPRC_r2_assemblies_0.6.1.agc \
  -n 233 \
  --subset-sequence-list subset.txt
```

The subset file must contain sequence names recognized by `impg similarity`.

## Progress and debug output

- default mode is silent
- `--progress` prints wrapper and `impg` progress messages
- `--debug` preserves pair-level intermediate files

Example:

```bash
bash scripts/ld/run_ld_pipeline.sh \
  --bed regions.bed \
  -p ../data/hprc465vschm13.aln.paf.gz \
  -f ../data/HPRC_r2_assemblies_0.6.1.agc \
  -n 233 \
  --progress \
  --debug
```

## Manual step-by-step workflow

The wrapper automates the following steps for a single pair of loci.

### 1. Similarity at locus 1

```bash
impg similarity \
  -r CHM13#0#chr22:11582486-11583486 \
  -a ../data/hprc465vschm13.aln.paf.gz \
  --sequence-files ../data/HPRC_r2_assemblies_0.6.1.agc \
  > loc1.sim
```

### 2. Similarity at locus 2

```bash
impg similarity \
  -r CHM13#0#chr22:11586486-11589486 \
  -a ../data/hprc465vschm13.aln.paf.gz \
  --sequence-files ../data/HPRC_r2_assemblies_0.6.1.agc \
  > loc2.sim
```

### 3. Joint haplotype matrix

```bash
python3 scripts/ld/cluster_two_loci.py \
  --locus1 loc1.sim \
  --locus2 loc2.sim \
  --haplotype-matrix-frequencies
```

This produces a normalized contingency table such as:

```text
loc1	l2c1	l2c2	l2c3	total
l1c1	0.80	0.10	0.00	0.90
l1c2	0.05	0.02	0.00	0.07
l1c3	0.01	0.01	0.01	0.03
total	0.86	0.13	0.01	1.00
```

### 4. Evaluate \(D_{ij}\) and summary statistics

```bash
python3 scripts/ld/compute_dij.py \
  --matrix loc1_loc2.haplotype_matrix.tsv \
  --output dij.tsv \
  --sample-size 233
```

The `dij.tsv` file contains the per-cell \(D_{ij}\) values, and the summary statistics are written to stderr by `compute_dij.py`.
