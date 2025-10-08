# How Hudson Fst is evaluated


EDAR 200kb , se possibbile allargare a 2Mb 
 echo -e "chr2\t109257703\t109457703" |  bedtools makewindows -b - -w 5000 > region.bed 

../../impop/scripts/run_h-fst.sh  -p ../../data/hprc465vschm13.aln.paf.gz -s ../../data/HPRC_r2_assemblies_0.6.1.agc  -b region.bed -A ../../metadata/agc.EUR   -B ../../metadata/agc.AFR  -o eur.afr.fst

../../impop/scripts/run_h-fst.sh  -p ../../data/hprc465vschm13.aln.paf.gz -s ../../data/HPRC_r2_assemblies_0.6.1.agc  -b region.bed -A ../../metadata/agc.EAS   -B ../../metadata/agc.AFR  -o eas.afr.fst

../../impop/scripts/run_h-fst.sh  -p ../../data/hprc465vschm13.aln.paf.gz -s ../../data/HPRC_r2_assemblies_0.6.1.agc  -b region.bed -A ../../metadata/agc.SAS   -B ../../metadata/agc.AFR  -o sas.afr.fst

../../impop/scripts/run_h-fst.sh  -p ../../data/hprc465vschm13.aln.paf.gz -s ../../data/HPRC_r2_assemblies_0.6.1.agc  -b region.bed -A ../../metadata/agc.AMR   -B ../../metadata/agc.AFR  -o amr.afr.fst

../../impop/scripts/run_h-fst.sh  -p ../../data/hprc465vschm13.aln.paf.gz -s ../../data/HPRC_r2_assemblies_0.6.1.agc  -b region.bed -A ../../metadata/agc.EAS   -B ../../metadata/agc.EUR  -o eas.eur.fst

../../impop/scripts/run_h-fst.sh  -p ../../data/hprc465vschm13.aln.paf.gz -s ../../data/HPRC_r2_assemblies_0.6.1.agc  -b region.bed -A ../../metadata/agc.SAS   -B ../../metadata/agc.EUR  -o sas.eur.fst


../../impop/scripts/run_h-fst.sh  -p ../../data/hprc465vschm13.aln.paf.gz -s ../../data/HPRC_r2_assemblies_0.6.1.agc  -b region.bed -A ../../metadata/agc.AMR   -B ../../metadata/agc.EUR  -o amr.eur.fst

../../impop/scripts/run_h-fst.sh  -p ../../data/hprc465vschm13.aln.paf.gz -s ../../data/HPRC_r2_assemblies_0.6.1.agc  -b region.bed -A ../../metadata/agc.EAS   -B ../../metadata/agc.SAS  -o eas.sas.fst

../../impop/scripts/run_h-fst.sh  -p ../../data/hprc465vschm13.aln.paf.gz -s ../../data/HPRC_r2_assemblies_0.6.1.agc  -b region.bed -A ../../metadata/agc.AMR   -B ../../metadata/agc.SAS  -o amr.sas.fst

../../impop/scripts/run_h-fst.sh  -p ../../data/hprc465vschm13.aln.paf.gz -s ../../data/HPRC_r2_assemblies_0.6.1.agc  -b region.bed -A ../../metadata/agc.AMR   -B ../../metadata/agc.EAS  -o amr.eas.fst


 ../../impop/scripts/plot_fst_trend.R --input-dir res  --title "EDAR " --dpi 300 -o edar.fst.png