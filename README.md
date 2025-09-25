software required [agc](https://github.com/refresh-bio/agc) impg, odgi 



### Nucleotide diversity 

Sample-level unbiased estimator of average pairwise nucleotide diversity, with corrections for finite sample size

\[
\hat{\pi} = \frac{n}{n-1} \sum_{ij} x_i x_j \pi_{ij} 
= \frac{n}{n-1} \sum_{i=2}^n \sum_{j=1}^{i-1} 2 x_i x_j \pi_{ij}
\]

Nei, M.; Li, W.-H. (1979). "Mathematical Model for Studying Genetic Variation in Terms of Restriction Endonucleases". *Proceedings of the National Academy of Sciences*. **76** (10): 5269–5273. doi:[10.1073/pnas.76.10.5269](https://doi.org/10.1073/pnas.76.10.5269). PMC [413122](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC413122). PMID [291943](https://pubmed.ncbi.nlm.nih.gov/291943/).


query the gfa to extract a window:
```
impg query -p ../data/hprc465vschm13.aln.paf.gz -r CHM13#0#chr1:158341439-158343639 --sequence-files ../data/HPRC_r2_assemblies_0.6.1.agc -o gfa >tmp.1.gfa && odgi sort -i tmp.1.gfa -o - | odgi view -i - -g >tmp.gfa && rm tmp.1.gfa
```

evaluate similarity across pairs of haplotypes in the window:
```
odgi similarity -i tmp.gfa > tmp.sim 
```

evaluate nucleotide diversity in the window 
```
python3 ../scr/pica2.2.py  tmp.sim  -t .988  -l 200 -r 5
```

the first two commands can be substituted bya single command [need a fix]: 
```
impg similarity -p hprc465vschm13.aln.paf.gz -r CHM13#0#chr1:158341439-158341639 --sequence-files HPRC_r2_assemblies_0.6.1.agc
```
