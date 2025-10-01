<img src="img/impop1.png" alt="implicit pangenome diagram" align="right" width="160" />

# Population Genomics Tools for Implicit Pangenomes

Software required: [agc](https://github.com/refresh-bio/agc), [impg](https://github.com/pangenome/impg), [odgi](https://github.com/pangenome/odgi), [povu](https://github.com/pangenome/povu)

[Dataset info](doc/where_hprc_data.md)


### Nucleotide diversity

Sample-level unbiased estimator of average pairwise nucleotide diversity, with corrections for finite sample size [wiki](https://en.wikipedia.org/wiki/Nucleotide_diversity):

$$\hat{\pi} = \frac{n}{n-1} \sum_{ij} x_i x_j \pi_{ij} = \frac{n}{n-1} \sum_{i=2}^n \sum_{j=1}^{i-1} 2 x_i x_j \pi_{ij}$$

> Nei, M.; Li, W.-H. (1979). "Mathematical Model for Studying Genetic Variation in Terms of Restriction Endonucleases". *Proceedings of the National Academy of Sciences*. **76** (10): 5269–5273. doi:[10.1073/pnas.76.10.5269](https://doi.org/10.1073/pnas.76.10.5269). PMC [413122](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC413122). PMID [291943](https://pubmed.ncbi.nlm.nih.gov/291943/).

#### [How to evaluate pi?](doc/how_pi.md)

### Hudson's Fst 

Measure of genetic differentiation between populations, comparing average diversity within and between them [wiki](https://en.wikipedia.org/wiki/Fixation_index):

$$F_{ST} = \frac{\pi_{between} - \pi_{within}}{\pi_{between}}$$

> Hudson, R. R.; Slatkin, M.; Maddison, W. P. (1992). "Estimation of levels of gene flow from DNA sequence data." *Genetics*. **132** (2): 583–589.


#### [How to evaluate Fst?](doc/how_tjd.md)

### Tajima's D 

Test statistic comparing two estimates of genetic variation: the average pairwise nucleotide diversity (π) and the number of segregating sites (S). It detects departures from the neutral mutation–drift equilibrium, such as population expansion, bottlenecks, or selection [wiki](https://en.wikipedia.org/wiki/Tajima%27s_D):

$$D = \frac{\pi - S/a_1}{\sqrt{e_1 S + e_2 S (S - 1)}}$$

**Interpretation**:  
- **D ≈ 0** → Consistent with neutral evolution and constant population size.  
- **D < 0** → Excess of low-frequency polymorphisms, suggesting population expansion or purifying selection.  
- **D > 0** → Excess of intermediate-frequency polymorphisms, suggesting balancing selection or a recent population bottleneck.  

> Tajima, F. (1989). *Statistical method for testing the neutral mutation hypothesis by DNA polymorphism.* Genetics, **123**(3), 585–595.  [Link to article (Genetics, 1989)](https://www.genetics.org/content/123/3/585)

##### Example

The Duffy antigen locus (*ACKR1* aka *DARC*) represents one of the most striking examples of natural selection in human evolution, driven by malaria resistance. Three major allelic variants exist: FYB (ancestral), FYA (common in Asia and Europe), and FYO (Duffy null, fixed in sub-Saharan Africa). The FYO allele, characterized by a promoter mutation that prevents DARC expression on red blood cells, provides near-complete protection against Plasmodium vivax malaria and has undergone one of the strongest selective sweeps in the human genome with a selection coefficient of 0.043. This ancient selective event began approximately 42,000 years ago from standing variation at very low frequency (0.1\%), rather than from a new mutation, and swept to near-fixation throughout equatorial Africa where P. vivax posed the greatest threat. The extreme geographic differentiation of these alleles - with FY\*O at >99% frequency in most sub-Saharan populations but virtually absent elsewhere - combined with signatures of reduced diversity and extended haplotype homozygosity, demonstrates how pathogen pressure has profoundly shaped human genetic variation. Interestingly, despite being a textbook example of positive selection, the complex evolutionary history of this locus, involving selection on standing variation rather than a simple hard sweep, means it is often missed by standard genome-wide selection scans.

> McManus, K. F.; et al. (2017). "Population genetic analysis of the DARC locus (Duffy) reveals adaptation from standing variation associated with malaria resistance in humans." *PLoS Genetics*. **13** (3): e1006560. doi:[10.1371/journal.pgen.1006560](https://doi.org/10.1371/journal.pgen.1006560).

> Voight, B. F.; Kudaravalli, S.; Wen, X.; Pritchard, J. K. (2006). "A map of recent positive selection in the human genome." *PLoS Biology*. **4** (3): e72. doi:[10.1371/journal.pbio.0040072](https://doi.org/10.1371/journal.pbio.0040072).


Evaluate nucleotide diversity in the ACKR1 gene region (chr1:158,340,000-158,344,000) in windows of 200 bp. Use the HPRCv2 assembly; coordinates are relative to CHM13.

Key SNP: rs2814778 (T-42C) defines the FY\*O (also known as Duffy null) in the GATA-1 transcription factor binding site in the *ACKR1* gene promoter region. 

| Interval | **hg38** | **CHM13** |
|----------|----------|-----------|
| SNP Position | chr1:159204893-159204893 (rs2814778) | chr1:158341919-158341919 (rs2814778) |
| 50 kb Interval | chr1:159179893-159229893 | chr1:158316925-158366917 |
| 100 kb Interval | chr1:159154893-159254893 | chr1:158291925-158391926 |
| 200 kb Interval | chr1:159104893-159304893 | chr1:158241938-158441927 |




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
