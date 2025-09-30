# How nucleotide diversity is evaluated

## impg similarity + pica2.py

##### One window

1. Generate the similarity matrix (requires impg support for AGC archives; adjust paths as needed):
```
impg similarity \
  -p hprc465vschm13.aln.paf.gz \
  -r CHM13#0#chr1:158341439-158341639 \
  --sequence-files HPRC_r2_assemblies_0.6.1.agc > tmp.sim
```

2. Evaluate nucleotide diversity for the window (recommended `-t 0.999`, `-r 4`):
```
python3 scripts/pica2.py tmp.sim \
  -t 0.999 
  -l 200 
  -r 5
```

##### One window with subsetting sequences

Generate the similarity matrix from a subset of assemblies:
```
impg similarity \
  -p hprc465vschm13.aln.paf.gz 
  -r CHM13#0#chr1:158341439-158341639 
  --sequence-files HPRC_r2_assemblies_0.6.1.agc 
  --subset-sequence-list /metadata/agc.EUR > EUR.sim
```

`agc.EUR` is a plain-text file (one assembly name per line)


##### Multi-window: run_pica2_impg.sh 

1. Prepare a BED file with windows (window size max 10kb as per impg similarity requirements):
```
echo -e "chr1\t158341439\t158341839" | bedtools makewindows -b - -w 200 > regions.bed
```

2. Run the wrapper (recommended `-t 0.999`, `-r 4`):
```
../impop/scripts/run_pica2_impg.sh -b regions.bed -t 0.999 -r 4
```

Use `-u` to restrict to assemblies listed in a plain-text file (one assembly name per line):
```
../impop/scripts/run_pica2_impg.sh -b regions.bed -t 0.999 -r 4 \
  -p hprc465vschm13.aln.paf.gz \
  -s HPRC_r2_assemblies_0.6.1.agc \
  -u /metadata/agc.EUR
```

Add `-l <length>` when you need to override the window length passed to `pica2.py` (defaults to the BED interval size).
Check `../impop/scripts/run_pica2_impg.sh -h` for the full option list.


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
Rscript scripts/plot_pi_trend.R \
  --input EUR=results/eur.pi.tsv \
  --input AFR=results/afr.pi.tsv \
  --title "ACKR1 nucleotide diversity" \
  --highlight chr1:158341439-158341639 \
  --output ackr1_pi_trend.png
```






## impg query + odgi similarity + pica2.py

###### one window 
1.query the gfa to extract a window:
```
impg query -p ../data/hprc465vschm13.aln.paf.gz -r CHM13#0#chr1:158341439-158343639 --sequence-files ../data/HPRC_r2_assemblies_0.6.1.agc -o gfa >tmp.1.gfa && odgi sort -i tmp.1.gfa -o - | odgi view -i - -g >tmp.gfa && rm tmp.1.gfa
```

2. evaluate similarity across pairs of haplotypes in the window:
```
odgi similarity -i tmp.gfa > tmp.sim 
```

3. evaluate nucleotide diversity in the window 
```
python3 scripts/pica2.py tmp.sim -t 0.988 -l 200 -r 5
```

###### run_pica2_odgi.sh: multi-window   
0. make a bed file for windows 
```
echo -e "chr1\t158340000\t158344000" | bedtools  makewindows -b - -w 200   > ackr1.win.bed
```

1. run the wrap (reccomended t 0.999, r 4 )
```
../impop/scripts/run_pica2_odgi.sh -b ackr1.win.bed  -t 0.999 -r 4
```