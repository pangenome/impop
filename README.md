population genomics tools for implicit pangenomes 


software required [agc](https://github.com/refresh-bio/agc) impg, odgi 



### Nucleotide diversity 

Sample-level unbiased estimator of average pairwise nucleotide diversity, with corrections for finite sample size [wiki](https://en.wikipedia.org/wiki/Nucleotide_diversity): 

$$\hat{\pi} = \frac{n}{n-1} \sum_{ij} x_i x_j \pi_{ij} = \frac{n}{n-1} \sum_{i=2}^n \sum_{j=1}^{i-1} 2 x_i x_j \pi_{ij}$$

>> Nei, M.; Li, W.-H. (1979). "Mathematical Model for Studying Genetic Variation in Terms of Restriction Endonucleases". *Proceedings of the National Academy of Sciences*. **76** (10): 5269–5273. doi:[10.1073/pnas.76.10.5269](https://doi.org/10.1073/pnas.76.10.5269). PMC [413122](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC413122). PMID [291943](https://pubmed.ncbi.nlm.nih.gov/291943/).


##### example 


Eavluate nucleotide diversity in the ACKR1 (DARC) gene region (chr1:158,340,000-158,344,000) in windows of 200bp
Use the HPRCv2, coordinates are relative to chm13 

### impg query + odgi similarity + pica2.py

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
python3 ../scr/pica2.2.py  tmp.sim  -t .988  -l 200 -r 5
```

###### wrap_pica2_impg_odgi.sh: multi-window   
0. make a bed file for windows 
```
echo -e "chr1\t158340000\t158344000" | bedtools  makewindows -b - -w 200   > ackr1.win.bed
```

1. run the wrap (reccomended t 0.999, r 4 )
```
../impop/scr/wrap_pica2_impg_odgi.sh -b ackr1.win.bed  -t 0.999 -r 4
```


### impg similarity + pica2.py

###### one window 
1-2.the first two commands can be substituted by a single command [need a fix for impg to work with agc]: 
```
impg similarity -p hprc465vschm13.aln.paf.gz -r CHM13#0#chr1:158341439-158341639 --sequence-files HPRC_r2_assemblies_0.6.1.agc
```

3. evaluate nucleotide diversity in the window 
```
python3 ../scr/pica2.2.py  tmp.sim  -t .988  -l 200 -r 5
```

###### wrap_pica2_impg.sh: multi-window  
0. make a bed file for windows 
```
echo -e "chr1\t158340000\t158344000" | bedtools  makewindows -b - -w 200   > ackr1.win.bed
```

 1. run the wrap (reccomended t 0.999, r 4 )
```
../impop/scr/wrap_pica2_impg.sh -b ackr1.win.bed  -t 0.999 -r 4
```


###### one window with subsetting sequences 
```
impg similarity -p ../data/hprc465vschm13.aln.paf.gz -r CHM13#0#chr1:158341439-158341639 --sequence-files ../data/HPRC_r2_assemblies_0.6.1.agc  --subset-sequence-list  ../metadata/agc.EUR  > EUR.sim 
```
