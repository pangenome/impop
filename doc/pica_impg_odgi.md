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