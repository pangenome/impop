# How nucleotide diversity is evaluated

## impg similarity + pica2.py
test region (200bp) chr1:158,341,439-158,341,639

##### One window

1. Generate the similarity matrix comparing all sequences:
```
impg similarity \
  -p hprc465vschm13.aln.paf.gz \
  -r CHM13#0#chr1:158341439-158341639 \
  --sequence-files HPRC_r2_assemblies_0.6.1.agc > tmp.sim
```

2. Evaluate nucleotide diversity for the window (recommended `-t 0.999`, `-r 5`):
```
python3 pica2.py tmp.sim \
  -t 0.999 
  -l 200 
  -r 5
```
>> pi for all sequences" 0.00000272 (sequence length: 200)


##### One window with subsetting sequences

Generate the similarity matrix from a subset of assemblies (`agc.EUR` is a plain-text file (one assembly name per line):
```
impg similarity \
  -p hprc465vschm13.aln.paf.gz 
  -r CHM13#0#chr1:158341439-158341639 
  --sequence-files HPRC_r2_assemblies_0.6.1.agc 
  --subset-sequence-list /metadata/agc.EUR > EUR.sim
```
>> 0.00000311 (sequence length: 200)


##### Multi-window: run_pica2_impg.sh 

1. Prepare a BED file with windows (window size max 10kb as per impg similarity requirements):
```
echo -e "chr1\t158341239\t158341839" | bedtools makewindows -b - -w 200 > regions.bed
```

2. Run the wrapper (recommended `-t 0.999`, `-r 4`):
Use `-u` to restrict to assemblies listed in a plain-text file (one assembly name per line)
Add `-l <length>` when you need to override the window length passed to `pica2.py` (defaults to the BED interval size).
Check `../impop/scripts/run_pica2_impg.sh -h` for the full option list.

```
run_pica2_impg.sh -b regions.bed -t 0.999 -r 4 \
  -p hprc465vschm13.aln.paf.gz \
  -s HPRC_r2_assemblies_0.6.1.agc \
  -u ../metadata/agc.EUR \
  -o pi.eur.tsv
```
>>REGION  SUBSET  LENGTH  THRESHOLD       R_VALUE PICA_OUTPUT
>>CHM13#0#chr1:158341239-158341439        agc.EUR 200     0.999   4       0.00000721 (sequence length: 200)
>>CHM13#0#chr1:158341439-158341639        agc.EUR 200     0.999   4       0.00000311 (sequence length: 200)
>>CHM13#0#chr1:158341639-158341839        agc.EUR 200     0.999   4       0.00000000 (sequence length: 200)


```
run_pica2_impg.sh -b regions.bed -t 0.999 -r 4 \
  -p hprc465vschm13.aln.paf.gz \
  -s HPRC_r2_assemblies_0.6.1.agc \
  -u ../metadata/agc.AFR \
  -o pi.afr.tsv
```
>>REGION  SUBSET  LENGTH  THRESHOLD       R_VALUE PICA_OUTPUT
>>CHM13#0#chr1:158341239-158341439        agc.AFR 200     0.999   4       0.00000639 (sequence length: 200)
>>CHM13#0#chr1:158341439-158341639        agc.AFR 200     0.999   4       0.00000488 (sequence length: 200)
>>CHM13#0#chr1:158341639-158341839        agc.AFR 200     0.999   4       0.00000000 (sequence length: 200)


### Plotting pi trends

Use `scripts/plot_pi_trend.R` to turn one or more `pica2.py` summary tables into a comparative trend plot.

Dependencies:
- R (tested with >= 4.2)
- R packages `ggplot2` and `dplyr`

The script expects tab-delimited inputs with at least `REGION` and `PICA_OUTPUT` columns (as emitted by `scripts/pica2.py`). Supply each table with `--input`, optionally prefixing a label (`--input EUR=results/eur.pi.tsv`). When no label is provided, the script derives one from the `SUBSET` column or the file name.

Common options:
- `--output` path for the plot (default `pi_trend.png`)
- `--title` custom plot title
- `--dpi` image resolution in dots per inch (default 150)
- `--highlight chrom:start-end` to shade specific intervals (repeatable)
- `--highlight-bed file.bed` to add intervals from a BED file

Example:
```
Rscript plot_pi_trend.R \
  --input EUR=pi.eur.tsv \
  --input AFR=pi.afr.tsv \
  --title "ACKR1" \
  --highlight chr1:158341439-158341639 \
  --output ackr1_pi.png
```

