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


#### [How to evaluate Fst?](doc/how_fst.md)

### Tajima's D 

Test statistic comparing two estimates of genetic variation: the average pairwise nucleotide diversity (π) and the number of segregating sites (S). It detects departures from the neutral mutation–drift equilibrium, such as population expansion, bottlenecks, or selection [wiki](https://en.wikipedia.org/wiki/Tajima%27s_D):

$$D = \frac{\pi - S/a_1}{\sqrt{e_1 S + e_2 S (S - 1)}}$$

**Interpretation**:  
- **D ≈ 0** → Consistent with neutral evolution and constant population size.  
- **D < 0** → Excess of low-frequency polymorphisms, suggesting population expansion or purifying selection.  
- **D > 0** → Excess of intermediate-frequency polymorphisms, suggesting balancing selection or a recent population bottleneck.  

> Tajima, F. (1989). *Statistical method for testing the neutral mutation hypothesis by DNA polymorphism.* Genetics, **123**(3), 585–595.  [Link to article (Genetics, 1989)](https://www.genetics.org/content/123/3/585)

#### [How to evaluate Fst?](doc/how_tjd.md)

### Example: positive selection at *EDAR* 

The *EDAR* (chr2:109355029-109449802 on CHM13) 370A variant shows one of the strongest signals of positive selection in the human genome. This mutation arose ~10,740 years ago and swept to near-fixation in East Asian and Native American populations while remaining virtually absent in Africans. The variant shows extreme population differentiation (Fst = 0.760) and enhances NF-κB activation, affecting the development of hair, teeth, and sweat glands.

The specific functional variant under selection is EDAR V370A (**rs3827760**, chr2:108897145-108897145 on hg38, chr2:109357703-109357703 on CHM13). This is a valine to alanine substitution at position 370 of the EDAR protein. The variant is located within the death domain of EDAR. 

| Interval | **hg38** | **CHM13** |
|----------|----------|-----------|
| SNP Position | chr2:108897145-108897145 | chr2:109357703-109357703 |
| 50 kb Interval | chr2:108872145-108922145 | chr2:109332703-109382703 |
| 100 kb Interval | chr2:108847145-108947145 | chr2:109307703-109407703 |
| 200 kb Interval | chr2:108797145-108997145 | chr2:109257703-109457703 |


> Bryk, J.; Hardouin, E.; Pugach, I.; Hughes, D.; Strotmann, R.; Stoneking, M.; Myles, S. (2008). "Positive Selection in East Asians for an EDAR Allele that Enhances NF-κB Activation". *PLOS ONE*. **3** (5): e2209. doi:[10.1371/journal.pone.0002209](https://doi.org/10.1371/journal.pone.0002209). PMC [2374902](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2374902). PMID [18493316](https://pubmed.ncbi.nlm.nih.gov/18493316/).

> Sabeti, P.C.; Varilly, P.; Fry, B.; Lohmueller, J.; Hostetter, E.; Cotsapas, C.; Xie, X.; Byrne, E.H.; McCarroll, S.A.; Gaudet, R.; Schaffner, S.F.; Lander, E.S.; The International HapMap Consortium (2007). "Genome-wide detection and characterization of positive selection in human populations". *Nature*. **449** (7164): 913–918. doi:[10.1038/nature06250](https://doi.org/10.1038/nature06250). PMC [2687721](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2687721). PMID [17943131](https://pubmed.ncbi.nlm.nih.gov/17943131/).


### Example: positive selection at *ACKR1* 

The Duffy antigen locus (ACKR1/DARC) exemplifies powerful natural selection through malaria resistance. The FYO allele, carrying a promoter mutation that blocks DARC expression on red blood cells, provides near-complete protection against Plasmodium vivax malaria. This variant underwent one of the strongest selective sweeps in humans (selection coefficient 0.043), rising from 0.1% frequency ~42,000 years ago to >99% in sub-Saharan Africa while remaining virtually absent elsewhere. Despite being a textbook case of positive selection, its evolution from standing variation rather than a new mutation causes it to be overlooked by standard selection scans.


Key SNP: rs2814778 (T-42C) defines the FY\*O (also known as Duffy null) in the GATA-1 transcription factor binding site in the *ACKR1* gene promoter region. 

| Interval | **hg38** | **CHM13** |
|----------|----------|-----------|
| SNP Position | chr1:159204893-159204893 (rs2814778) | chr1:158341919-158341919 (rs2814778) |
| 50 kb Interval | chr1:159179893-159229893 | chr1:158316925-158366917 |
| 100 kb Interval | chr1:159154893-159254893 | chr1:158291925-158391926 |
| 200 kb Interval | chr1:159104893-159304893 | chr1:158241938-158441927 |


> McManus, K. F.; et al. (2017). "Population genetic analysis of the DARC locus (Duffy) reveals adaptation from standing variation associated with malaria resistance in humans." *PLoS Genetics*. **13** (3): e1006560. doi:[10.1371/journal.pgen.1006560](https://doi.org/10.1371/journal.pgen.1006560).

> Voight, B. F.; Kudaravalli, S.; Wen, X.; Pritchard, J. K. (2006). "A map of recent positive selection in the human genome." *PLoS Biology*. **4** (3): e72. doi:[10.1371/journal.pbio.0040072](https://doi.org/10.1371/journal.pbio.0040072).

