# How Fst is evaluated

##### One window

1. Generate a 500 bp window:
```
echo -e "chr2\t109332703\t109382703" | bedtools makewindows -b - -w 500 > darc.bed
```

2. Evaluate nucleotide diversity for subset A (`piA`):
```
impg similarity \
  --sequence-files ../../data/HPRC_r2_assemblies_0.6.1.agc \
  -p ../../data/hprc465vschm13.aln.paf.gz \
  -r CHM13#0#chr1:158291925-158292425 \
  --subset-sequence-list ../../metadata/agc.EUR > eur.sim

python3 ../../impop/scripts/pica2.py eur.sim \
  -t 0.999 \
  -l 500 \
  -r 5
# -> 0.00510773 (sequence length: 500)
```

3. Evaluate nucleotide diversity for subset B (`piB`):
```
impg similarity \
  --sequence-files ../../data/HPRC_r2_assemblies_0.6.1.agc \
  -p ../../data/hprc465vschm13.aln.paf.gz \
  -r CHM13#0#chr1:158291925-158292425 \
  --subset-sequence-list ../../metadata/agc.AFR > afr.sim

python3 ../../impop/scripts/pica2.py afr.sim \
  -t 0.999 \
  -l 500 \
  -r 5
# -> 0.15271961 (sequence length: 500)
```

4. Build the union list and evaluate the combined diversity (`piC`):
```
cat ../../metadata/agc.AFR ../../metadata/agc.EUR > temp.agc.AFREUR

impg similarity \
  --sequence-files ../../data/HPRC_r2_assemblies_0.6.1.agc \
  -p ../../data/hprc465vschm13.aln.paf.gz \
  -r CHM13#0#chr1:158291925-158292425 \
  --subset-sequence-list temp.agc.AFREUR > afreur.sim

python3 ../../impop/scripts/pica2.py eurafr.sim \
  -t 0.999 \
  -l 500 \
  -r 5
# -> 0.31322113 (sequence length: 500)

rm temp.agc.AFREUR
```

5. Compute the summary statistics:
```
piAB = 0.5 * (0.00510773 + 0.15271961)
Fst  = (0.31322113 - piAB) / 0.31322113
```




## run_fst_impg.sh

1. Prepare a BED file with windows that satisfy `impg similarity` constraints (recommended ≤10 kb):
```
echo -e "chr2\t109332703\t109382703" | bedtools makewindows -b - -w 5000 > region.bed
```

2. Ensure you have two subset lists (plain text, one assembly per line) describing unique populations A and B:
```
head ../../metadata/agc.EUR
head ../../metadata/agc.AFR
```

3. Run the wrapper, providing both subset lists, the BED file, and the `pica2.py` options (window length inferred from the BED entry):
```
../impop/scripts/run_fst_impg.sh \
  -A ../../metadata/agc.EUR \
  -B ../../metadata/agc.AFR \
  -b darc.bed \
  -p ../../data/hprc465vschm13.aln.paf.gz \
  -s ../../data/HPRC_r2_assemblies_0.6.1.agc \
  -t 0.999 \
  -r 5 \
  -o eur.afr.fst
```

4. Inspect the output table (`REGION`, `LENGTH`, `THRESHOLD`, `R_VALUE`, `PI_A`, `PI_B`, `PI_C`, `PI_AB_AVG`, `FST`). The script reports `NA` for `FST` when `piC` is zero and writes detailed `pica2.py` logs to `pica2_logs/` (override with `-d`).

Use `-P` to provide a different region prefix, and run `scripts/run_fst_impg.sh -h` for the full option list.

### Plotting Fst trends

Use `scripts/plot_fst_trend.R` to visualise windowed Fst estimates produced by `scripts/run_fst_impg.sh`.

Dependencies:
- R (tested with ≥ 4.2)
- R packages `ggplot2` and `dplyr`

The script expects tab-delimited inputs with `REGION` and `FST` columns. Provide each table with `--input`, optionally giving it a label (`--input AFRvEUR=results/darc.fst.tsv`). When no label is supplied the script derives one from the file name.

Common options mirror `plot_pi_trend.R`:
- `--output` path for the plot (default `fst_trend.png`)
- `--title` custom plot title
- `--dpi` image resolution (default 150)
- `--highlight chrom:start-end` to shade specific intervals (repeatable)
- `--highlight-bed file.bed` to add intervals from a BED file

Example:
```
Rscript scripts/plot_fst_trend.R \
  --input AFRvEUR=results/darc.fst.tsv \
  --title "ACKR1 Fst" \
  --highlight chr1:158291925-158292425 \
  --output ackr1_fst_trend.png
```
