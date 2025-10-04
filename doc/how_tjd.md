# How Tajima's D is evaluated


###### one window 

1. choose a window: chr1:158341439-158343639 

2. determine the number of segregating sites (S) 
2a. query the gfa to extract the window :
```
impg query -p ../data/hprc465vschm13.aln.paf.gz -r CHM13#0#chr1:158341439-158343639 --sequence-files ../data/HPRC_r2_assemblies_0.6.1.agc -o gfa >tmp.1.gfa && odgi sort -i tmp.1.gfa -o - | odgi view -i - -g >tmp.gfa && rm tmp.1.gfa
```
2b. use the gfa to evaluate number of segregating sites (S)
```
povu gfa2vcf  -i tmp.gfa  --stdout  CHM13 | grep -v "^#" | wc -l
```
>>> S is the standard out of povu 

3. determine nucleotide diversity (pi) 
3a. Generate the similarity matrix (requires impg support for AGC archives; adjust paths as needed):
```
impg similarity \
  -p hprc465vschm13.aln.paf.gz \
  -r CHM13#0#chr1:158341439-158341639 \
  --sequence-files HPRC_r2_assemblies_0.6.1.agc > tmp.sim
```

3b. Evaluate nucleotide diversity for the window (recommended `-t 0.999`, `-r 4`):
```
python3 ../impop/scripts/pica2.py tmp.sim \
  -t 0.999 
  -l 200 
  -r 5
```
>>> pi  is  the result of pica2.py 

4. evaluate number of sequences  (n)
```
wc -l all.agc
```
>>> n is the number of lines in the sample list file  

4. Evaluate Tajimas'D  using n, p ans S 

python3 ../impop/scripts/tj_d.py  -n 446 -p 0.59146123 -S 20


##### Multi-window: run_pica2_impg.sh

1. Prepare a BED file of windows (≤10 kb for `impg similarity`):
```
echo -e "chr2\t109332703\t109382703" | bedtools makewindows -b - -w 5000 > region.bed
```

2. 

```
../impop/scripts/run_tajd.sh \
  -b regions.bed \
  -l ../metadata/all.agc \
  -p ../data/hprc465vschm13.aln.paf.gz \
  -s ../data/HPRC_r2_assemblies_0.6.1.agc \
  -t 0.999 \
  -r 5 \
  -o regions.tajd.tsv
```

The output table reports, per window, the observed segregating sites (S), mean pairwise diversity (π), sample count (n), and Tajima's D.




### Tajima's D per window

Use `scripts/run_tajd.sh` to compute segregating sites (S), nucleotide diversity (π), sample count (n), and Tajima's D for each BED window by combining `impg query`, `odgi`, `povu gfa2vcf`, `impg similarity`, `scripts/pica2.py`, and `scripts/tj_d.py`.

Required inputs:
- `-b` BED file with windows
- `-l` sample list (one sequence ID per line) used as the subset for all analyses

Common options:
- `-p` PAF alignment (`impg query/similarity`)
- `-s` AGC archive of assemblies
- `-t` / `-r` parameters forwarded to `scripts/pica2.py`
- `-P` region prefix (default `CHM13#0#`) and `-R` reference name passed to `povu`
- `-o` output TSV path (defaults to stdout)

Example:
```
scripts/run_tajd.sh \
  -b darc.bed \
  -l ../../metadata/all.agc \
  -p ../../data/hprc465vschm13.aln.paf.gz \
  -s ../../data/HPRC_r2_assemblies_0.6.1.agc \
  -t 0.999 \
  -r 5 \
  -o darc.tajd.tsv
```

The resulting table reports `REGION`, window `LENGTH`, number of `SAMPLES`, segregating sites (`SEGREGATING_SITES`), window-wide `PI`, and `TAJIMAS_D` (with zero-S windows yielding `NA`).

### Plotting trends across runs

Use `scripts/plot_tajd_trend.R` to visualise Tajima's D profiles from one or more `run_tajd.sh` outputs. Supply each file with `--input`, optionally prefixing a label before the equals sign. You can also highlight genomic intervals via `--highlight chrom:start-end` or `--highlight-bed path/to/regions.bed`.

Example:
```
Rscript scripts/plot_tajd_trend.R \
  --input=Baseline=results/baseline.tajd.tsv \
  --input=Treatment=results/treatment.tajd.tsv \
  --output tajd_trend.png \
  --title "Tajima's D comparison" \
  --highlight chr2:109332703-109382703
```

The script combines windows from all runs, lays them out along a concatenated genome axis, and draws per-run lines and points for Tajima's D while shading any highlighted regions.
