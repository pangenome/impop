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
echo -e "chr1\t158341439\t158341839" | bedtools makewindows -b - -w 200 > regions.bed
```

2. Run the wrapper to generate π per window and intermediate similarity files:
```
../impop/scripts/run_pica2_impg.sh \
  -b regions.bed \
  -t 0.999 \
  -r 5 \
  -p ../data/hprc465vschm13.aln.paf.gz \
  -s ../data/HPRC_r2_assemblies_0.6.1.agc \
  -u ../metadata/all.agc
```

3. Feed the resulting π values, sample list, and segregating-site counts into `scripts/run_tajd.sh` to obtain Tajima's D across all windows:
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
