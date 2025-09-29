population genomics tools for implicit pangenomes

Software required:
- [agc](https://github.com/refresh-bio/agc)
- impg
- odgi

### Nucleotide diversity

Sample-level unbiased estimator of average pairwise nucleotide diversity, with corrections for finite sample size [wiki](https://en.wikipedia.org/wiki/Nucleotide_diversity):

$$\hat{\pi} = \frac{n}{n-1} \sum_{ij} x_i x_j \pi_{ij} = \frac{n}{n-1} \sum_{i=2}^n \sum_{j=1}^{i-1} 2 x_i x_j \pi_{ij}$$

>> Nei, M.; Li, W.-H. (1979). "Mathematical Model for Studying Genetic Variation in Terms of Restriction Endonucleases". *Proceedings of the National Academy of Sciences*. **76** (10): 5269–5273. doi:[10.1073/pnas.76.10.5269](https://doi.org/10.1073/pnas.76.10.5269). PMC [413122](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC413122). PMID [291943](https://pubmed.ncbi.nlm.nih.gov/291943/).

##### Example

Evaluate nucleotide diversity in the ACKR1 (DARC) gene region (chr1:158,340,000-158,344,000) in windows of 200 bp. Use the HPRCv2 assembly; coordinates are relative to CHM13.

### impg similarity + pica2.py

#### One window

1. Generate the similarity matrix (requires impg support for AGC archives; adjust paths as needed):
```
impg similarity -p hprc465vschm13.aln.paf.gz -r CHM13#0#chr1:158341439-158341639 --sequence-files HPRC_r2_assemblies_0.6.1.agc > tmp.sim
```

2. Evaluate nucleotide diversity for the window (adjust the path to `scripts/pica2.py` if you run the command from elsewhere):
```
python3 scripts/pica2.py tmp.sim -t 0.988 -l 200 -r 5
```


#### One window with subsetting sequences

Generate the similarity matrix from a subset of assemblies:
```
impg similarity -p ../data/hprc465vschm13.aln.paf.gz -r CHM13#0#chr1:158341439-158341639 --sequence-files ../data/HPRC_r2_assemblies_0.6.1.agc --subset-sequence-list ../metadata/agc.EUR > EUR.sim
```


#### run_pica2_impg.sh: multi-window

1. Prepare a BED file with windows (window size max 10kb as per impg similarity requirements):
```
echo -e "chr1\t158341439\t158341839" | bedtools makewindows -b - -w 200 > ackr1.win.bed
```

2. Run the wrapper (recommended `-t 0.999`, `-r 4`):
```
../impop/scripts/run_pica2_impg.sh -b ackr1.win.bed -t 0.999 -r 4
```

Use `-p` and `-s` to override the default PAF and sequence archives, and `-u` to restrict to assemblies listed in a plain-text file (one assembly name per line):
```
../impop/scripts/run_pica2_impg.sh -b ackr1.win.bed -t 0.999 -r 4 \
  -p ../data/hprc465vschm13.aln.paf.gz \
  -s ../data/HPRC_r2_assemblies_0.6.1.agc \
  -u ../metadata/agc.EUR
```

Add `-l <length>` when you need to override the window length passed to `pica2.py` (defaults to the BED interval size).

Check `../impop/scripts/run_pica2_impg.sh -h` for the full option list.

