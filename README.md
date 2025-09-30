<img src="img/impop1.png" alt="implicit pangenome diagram" align="right" width="160" />

# Population Genomics Tools for Implicit Pangenomes

Software required:
- [agc](https://github.com/refresh-bio/agc)
- [impg](https://github.com/pangenome/impg)
- [odgi](https://github.com/pangenome/odgi)
- [povu](https://github.com/pangenome/povu)


[Dataset info](doc/where_hprc_data.md)


### Nucleotide diversity

Sample-level unbiased estimator of average pairwise nucleotide diversity, with corrections for finite sample size [wiki](https://en.wikipedia.org/wiki/Nucleotide_diversity):

$$\hat{\pi} = \frac{n}{n-1} \sum_{ij} x_i x_j \pi_{ij} = \frac{n}{n-1} \sum_{i=2}^n \sum_{j=1}^{i-1} 2 x_i x_j \pi_{ij}$$

>> Nei, M.; Li, W.-H. (1979). "Mathematical Model for Studying Genetic Variation in Terms of Restriction Endonucleases". *Proceedings of the National Academy of Sciences*. **76** (10): 5269–5273. doi:[10.1073/pnas.76.10.5269](https://doi.org/10.1073/pnas.76.10.5269). PMC [413122](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC413122). PMID [291943](https://pubmed.ncbi.nlm.nih.gov/291943/).

#### [How was pi calculated?](doc/how_pi.md)

### Hudson's Fst 

Hudson, Richard R., Montgomery Slatkin, and Wayne P. Maddison. "Estimation of levels of gene flow from DNA sequence data." Genetics 132.2 (1992): 583-589.

Voight, Benjamin F., et al. "A map of recent positive selection in the human genome." PLoS biology 4.3 (2006): e72.

#### [How was Fst calculated?](doc/how_fst.md)

##### Example

The Duffy antigen locus (DARC) represents one of the most striking examples of natural selection in human evolution, driven by malaria resistance. Three major allelic variants exist: FYB (ancestral), FYA (common in Asia and Europe), and FYO (Duffy null, fixed in sub-Saharan Africa). The FYO allele, characterized by a promoter mutation that prevents DARC expression on red blood cells, provides near-complete protection against Plasmodium vivax malaria and has undergone one of the strongest selective sweeps in the human genome with a selection coefficient of 0.043. This ancient selective event began approximately 42,000 years ago from standing variation at very low frequency (0.1\%), rather than from a new mutation, and swept to near-fixation throughout equatorial Africa where P. vivax posed the greatest threat. The extreme geographic differentiation of these alleles - with FY\*O at >99% frequency in most sub-Saharan populations but virtually absent elsewhere - combined with signatures of reduced diversity and extended haplotype homozygosity, demonstrates how pathogen pressure has profoundly shaped human genetic variation. Interestingly, despite being a textbook example of positive selection, the complex evolutionary history of this locus, involving selection on standing variation rather than a simple hard sweep, means it is often missed by standard genome-wide selection scans.



Evaluate nucleotide diversity in the ACKR1 gene region (chr1:158,340,000-158,344,000) in windows of 200 bp. Use the HPRCv2 assembly; coordinates are relative to CHM13.

Key SNP: rs2814778 located within this interval

This spans the ACKR1 coding sequence plus regulatory regions implicated in the selective sweep.


FY\*O (also known as Duffy null) is defined by a mutation (T-42C, rs2814778) in the GATA-1 transcription factor binding site in the DARC gene promoter region,


**hg38** Human: Common (1000 Genomes Phase 3 MAF >= 1%) Short Genetic Variants from dbSNP Release 155 (rs2814778)
dbSNP: rs2814778 Position: chr1:159204893-159204893
chr1:159,179,893-159,229,893 (50kb interval surrounding rs2814778)
chr1:159,154,893-159,254,893 (100kb interval surrounding rs2814778)
chr1:159,104,893-159,304,893 (200kb interval surrounding rs2814778)

**CHM13**
chr1:158341919-158341919 rs2814778
chr1:158316925-158366917 (50kb interval surrounding rs2814778)
chr1:158291925-158391926 (100kb interval surrounding rs2814778)
chr1:158241938-158441927 (200kb interval surrounding rs2814778)

| Reference | SNP Position | 50kb Interval | 100kb Interval | 200kb Interval |
|-----------|--------------|---------------|----------------|----------------|
| **hg38** | chr1:159204893-159204893 (rs2814778) | chr1:159179893-159229893 | chr1:159154893-159254893 | chr1:159104893-159304893 |
| **CHM13** | chr1:158341919-158341919 (rs2814778) | chr1:158316925-158366917 | chr1:158291925-158391926 | chr1:158241938-158441927 |Retry



McManus, Kimberly F., et al. "Population genetic analysis of the DARC locus (Duffy) reveals adaptation from standing variation associated with malaria resistance in humans." PLoS genetics 13.3 (2017): e1006560.

Voight, Benjamin F., et al. "A map of recent positive selection in the human genome." PLoS biology 4.3 (2006): e72.


### impg similarity + pica2.py

##### One window

1. Generate the similarity matrix (requires impg support for AGC archives; adjust paths as needed):
```
impg similarity -p hprc465vschm13.aln.paf.gz -r CHM13#0#chr1:158341439-158341639 --sequence-files HPRC_r2_assemblies_0.6.1.agc > tmp.sim
```

2. Evaluate nucleotide diversity for the window (adjust the path to `scripts/pica2.py` if you run the command from elsewhere):
```
python3 scripts/pica2.py tmp.sim -t 0.988 -l 200 -r 5
```

##### One window with subsetting sequences

Generate the similarity matrix from a subset of assemblies:
```
impg similarity -p ../data/hprc465vschm13.aln.paf.gz -r CHM13#0#chr1:158341439-158341639 --sequence-files ../data/HPRC_r2_assemblies_0.6.1.agc --subset-sequence-list ../metadata/agc.EUR > EUR.sim
```

##### Multi-window: run_pica2_impg.sh 

1. Prepare a BED file with windows (window size max 10kb as per impg similarity requirements):
```
echo -e "chr1\t158341439\t158341839" | bedtools makewindows -b - -w 200 > ackr1.win.bed
```

2. Run the wrapper (recommended `-t 0.999`, `-r 4`):
```
scripts/run_pica2_impg.sh -b ackr1.win.bed -t 0.999 -r 4
```

Use `-p` and `-s` to override the default PAF and sequence archives, and `-u` to restrict to assemblies listed in a plain-text file (one assembly name per line):
```
scripts/run_pica2_impg.sh -b ackr1.win.bed -t 0.999 -r 4 \
  -p ../data/hprc465vschm13.aln.paf.gz \
  -s ../data/HPRC_r2_assemblies_0.6.1.agc \
  -u ../metadata/agc.EUR
```

Add `-l <length>` when you need to override the window length passed to `pica2.py` (defaults to the BED interval size).

Check `scripts/run_pica2_impg.sh -h` for the full option list.

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

### Estimating Fst from nucleotide diversity

`scripts/run_fst_impg.sh` reproduces the manual workflow of computing nucleotide diversity in two populations (A, B), their union, and deriving windowed Fst.

Steps performed per window:
1. Evaluate `pi_A` using subset list A.
2. Evaluate `pi_B` using subset list B.
3. Combine both lists, evaluate `pi_union`, and report `pi_mean = 0.5*(pi_A + pi_B)`.
4. Emit `Fst = (pi_union - pi_mean) / pi_union`.

Required arguments:
- `-b` BED file with genomic windows
- `-t` similarity threshold and `-r` r-parameter for `pica2.py`
- `-A` / `-B` disjoint subset lists (blank lines and comments are ignored)

Optional arguments:
- `-p` PAF path and `-s` AGC archive overrides (passed to `run_pica2_impg.sh`)
- `-l` explicit sequence length override
- `-o` output TSV path (defaults to stdout)
- `-k` keep the temporary working directory for debugging

Example matching the ACKR1/DARC region:

```
scripts/run_fst_impg.sh \
  -b darc.bed \
  -t 0.999 \
  -r 5 \
  -A ../../metadata/agc.EUR \
  -B ../../metadata/agc.AFR \
  -p ../../data/hprc465vschm13.aln.paf.gz \
  -s ../../data/HPRC_r2_assemblies_0.6.1.agc \
  -o darc.fst
```
