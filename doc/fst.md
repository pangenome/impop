Hudson, Richard R., Montgomery Slatkin, and Wayne P. Maddison. "Estimation of levels of gene flow from DNA sequence data." Genetics 132.2 (1992): 583-589.


https://www.biorxiv.org/content/10.1101/2024.09.24.614506v1

I want to evaluate fst within genomic regions 

a) establish the genomic region 
b) evaluate pi in A --> piA
c) evaluate pi in B  --> piB 
d) merge list A and B in a new list C that contains all elements of A and B 
e) evaluate pi in list C --> piC 
f) evaluate average of piA and piB --> 0.5*(piA +piB)
g) evaluate Fst as (piC-piAB)/piC 


>> this translate into this series of commands: 

a) echo -e "chr1\t158291925\t158292425" | bedtools makewindows -b - -w 500 > darc.bed

b) impg similarity  --sequence-files ../../data/HPRC_r2_assemblies_0.6.1.agc -p  ../../data/hprc465vschm13.aln.paf.gz -r CHM13#0#chr1:158291925-158292425 --subset-sequence-list ../../metadata/agc.EUR > eur.sim 
python3 ../../impop/scripts/pica2.py eur.sim -t  0.999 -l 500  -r 5 
0.00510773 (sequence length: 500)

c) impg similarity  --sequence-files ../../data/HPRC_r2_assemblies_0.6.1.agc -p  ../../data/hprc465vschm13.aln.paf.gz -r CHM13#0#chr1:158291925-158292425 --subset-sequence-list ../../metadata/agc.AFR > afr.sim 
python3 ../../impop/scripts/pica2.py afr.sim -t  0.999 -l 500  -r 5 
0.15271961 (sequence length: 500)

d) cat ../../metadata/agc.AFR ../../metadata/agc.EUR > temp.agc.AFREUR  

e) impg similarity  --sequence-files ../../data/HPRC_r2_assemblies_0.6.1.agc -p  ../../data/hprc465vschm13.aln.paf.gz -r CHM13#0#chr1:158291925-158292425 --subset-sequence-list temp.agc.AFREUR > afreur.sim  
python3 ../../impop/scripts/pica2.py eurafr.sim -t  0.999 -l 500  -r 5 
0.31322113 (sequence length: 500)
rm temp.agc.AFREUR

f) average_piAB = 0.5*(0.00510773+0.15271961)
g) Fst is then evaluated as (0.31322113 - average_piAB )/0.31322113 


inspired by run_pica2_impg.sh write a new wrapper that takes in input two subset_list list A and list B (each with a unique set of sequences), a bed file with genomic regions, a paf file, a sequence file, -t -r for pica2.py (interval length inferred from bed) and reproduce the above steps  