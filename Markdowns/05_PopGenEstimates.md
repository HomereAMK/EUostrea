PopGenestimates:
Thetas, Tajima, Neutrality tests estimations from angsd for each populations.
================

- [Thetas, Tajima, Neutrality tests estimations from angsd for each populations.](#thetas-tajima-neutrality-tests-estimations-from-angsd-for-each-populations)
  - [Load module angsd.](#load-module-angsd)
  - [Get lists of samples](#get-lists-of-samples)
  - [Runs ANGSD under -doSaf on all populations, estimate the site allele frequency likelihood for each pop.](#runs-angsd-under--dosaf-on-all-populations-estimate-the-site-allele-frequency-likelihood-for-each-pop)
  - [Runs realSFS](#runs-realsfs)
  - [Calculates the thetas for each site:](#calculates-the-thetas-for-each-site)
  - [Gets summary of thetas:](#gets-summary-of-thetas)
  - [Performs final calculations:](#performs-final-calculations)
  - [Gets thetas per sliding window:](#gets-thetas-per-sliding-window)
  - [Edits thetas per sliding window:](#edits-thetas-per-sliding-window)
  - [Heterozygosity:](#Heterozygosity)




## Load module angsd.
```
#angsd
module load tools computerome_utils/2.0
module load htslib/1.16
module load bedtools/2.30.0
module load pigz/2.3.4
module load parallel/20220422
module load angsd/0.937

#R
module load tools computerome_utils/2.0
module load perl/5.30.2
module load pigz/2.3.4
module load pandas-profiling/3.0.0
module load graph-tool/2.43
module load gcc/11.1.0
module load intel/perflibs/64/2020_update2
module load R/4.2.0
```

## Get lists of samples
```
POP=("MOLU" "ZECE" "CRES" "ORIS" "CORS" "PONT" "RIAE" "MORL" "USAM" "TOLL" "COLN" "BARR" "TRAL" "CLEW" "RYAN" "GREV" "WADD" "NISS" "LOGS" "VENO" "HALS" "THIS" "KALV" "HYPP" "LANG" "BUNN" "DOLV" "HAUG" "HAFR" "INNE" "VAGS" "AGAB" "OSTR")
BAMLIST=/home/projects/dp_00007/people/hmon/EUostrea/01_infofiles/bamlist_EUostrea.txt

for query in ${POP[*]}
do 
    grep ${query} $BAMLIST > /home/projects/dp_00007/people/hmon/EUostrea/01_infofiles/EUostrea_${query}-Fst.list

done
```
🤝

## Runs ANGSD under -doSaf on all populations, estimate the site allele frequency likelihood for each pop.
```
for query in ${POP[*]}
    do
        REF=/home/projects/dp_00007/people/hmon/AngsdPopStruct/01_infofiles/fileOegenome10scaffoldC3G.fasta
        /home/projects/dp_00007/apps/Scripts/wrapper_angsd.sh -debug 2 -nThreads 40 -ref $REF -anc $REF -bam /home/projects/dp_00007/people/hmon/EUostrea/01_infofiles/EUostrea_${query}-Fst.list -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 -minMapQ 20 -minQ 20 -setMinDepthInd 1 -setMinDepth 600 -setMaxDepth 120 -GL 1 -doSaf 1 -out /home/projects/dp_00007/data/hmon/angsd_PopGen/Jan23--Unfolded_PopGen_${query}
    done
done
```
for query in ${POP[*]}
    do
        REF=/home/projects/dp_00007/people/hmon/AngsdPopStruct/01_infofiles/fileOegenome10scaffoldC3G.fasta
        /home/projects/dp_00007/apps/Scripts/wrapper_angsd.sh -debug 2 -nThreads 40 -ref $REF -anc $REF -bam /home/projects/dp_00007/people/hmon/EUostrea/01_infofiles/EUostrea_${query}-Fst.list -remove_bads 1 -uniqueOnly 1 -baq 1 -C 50 -minMapQ 20 -minQ 20 -minInd $((N_IND*2/3)) -GL 1 -doSaf 1  -out /home/projects/dp_00007/data/hmon/angsd_PopGen/Jan23--Unfolded_PopGenGEO_${query}
         done
done

## Runs realSFS
```
POP=("MOLU" "ZECE" "CRES" "ORIS" "CORS" "PONT" "RIAE" "MORL" "USAM" "TOLL" "COLN" "BARR" "TRAL" "CLEW" "RYAN" "GREV" "WADD" "NISS" "LOGS" "VENO" "HALS" "THIS" "KALV" "HYPP" "LANG" "BUNN" "DOLV" "HAUG" "HAFR" "INNE" "VAGS" "AGAB" "OSTR")

for query in ${POP[*]}

do
    realSFS -P 40 -fold 1 /home/projects/dp_00007/data/hmon/angsd_PopGen/Dec22--Unfolded_PopGen_${query}.saf.idx > /home/projects/dp_00007/data/hmon/angsd_PopGen/Dec22--Ind581_${query}.sfs
    
done
```

## Calculates the thetas for each site:
```
for query in ${POP[*]}

do
    realSFS saf2theta /home/projects/dp_00007/data/hmon/angsd_PopGen/Dec22--Unfolded_PopGen_${query}.saf.idx -sfs /home/projects/dp_00007/data/hmon/angsd_PopGen/Dec22--Ind581_${query}.sfs -fold 1 -outname /home/projects/dp_00007/data/hmon/angsd_PopGen/Dec22--Ind581_${query}_PopGenEstimates

done
```


## Gets summary of thetas:
```
for query in ${POP[*]}

do
    thetaStat print /home/projects/dp_00007/data/hmon/angsd_PopGen/Dec22--Ind581_${query}_PopGenEstimates.thetas.idx > /home/projects/dp_00007/data/hmon/angsd_PopGen/Dec22--Ind581_${query}_PopGenEstimates.Print
done
```

## Performs final calculations:
```
for query in ${POP[*]}

do
    N_IND=`cat /home/projects/dp_00007/people/hmon/EUostrea/01_infofiles/EUostrea_${query}-Fst.list | wc -l`
    Rscript --vanilla --slave /home/projects/dp_00007/apps/Scripts/GetsThetaSummaries.R /home/projects/dp_00007/data/hmon/angsd_PopGen/Dec22--Ind581_${query}_PopGenEstimates.Print $N_IND $query
done > /home/projects/dp_00007/data/hmon/angsd_PopGen/Dec22--Ind581.PopGenEstimates.txt
```

## Gets thetas per sliding window:
```
for query in ${POP[*]}

do
    thetaStat do_stat /home/projects/dp_00007/data/hmon/angsd_PopGen/Dec22--Ind581_${query}_PopGenEstimates.thetas.idx -win 20000 -step 20000 -outnames /home/projects/dp_00007/data/hmon/angsd_PopGen/Dec22--Ind581__20K_${query}_PopGenEstimates-Windows
done
```

## Edits thetas per sliding window:
```
for query in ${POP[*]}

do
    cut -f 2,3,4,5,9,14 /home/projects/dp_00007/data/hmon/angsd_PopGen/Dec22--Ind581__20K_${query}_PopGenEstimates-Windows.pestPG | tail -n +2 | sed -r 's/LG//g' | sed 's/^0*//' | awk '$6 > 0' | awk '{print $1"\t"$1":"$2"\t"$2-20000"\t"$2"\t"$6"\t"$3"\t"$3/$6"\t"$4"\t"$4/$6"\t"$5}' | awk 'BEGIN{print "CHR\tSNP\tgPoint\tWindow\tNumberOfSites\tsumTw\tTw\tsumTp\tTp\tTd"}1' > /home/projects/dp_00007/data/hmon/angsd_PopGen/Dec22--Ind581__20K_${query}_PopGenEstimates-Windows.tsv
done
```


## Heterozygosity
# Runs _ANGSD
```
N_IND=`cat /home/projects/dp_00007/people/hmon/EUostrea/01_infofiles/bamlist_EUostrea.txt | wc -l`
REF=/home/projects/dp_00007/people/hmon/AngsdPopStruct/01_infofiles/fileOegenome10scaffoldC3G.fasta
BAMLIST=/home/projects/dp_00007/people/hmon/EUostrea/01_infofiles/bamlist_EUostrea.txt



/home/projects/dp_00007/apps/Scripts/wrapper_angsd.sh \
-debug 2 -nThreads 40 -ref $REF -bam $BAMLIST \

/home/projects/dp_00007/apps/Scripts/wrapper_angsd.sh -debug 2 -nThreads 40 -ref $REF -bam $BAMLIST -remove_bads 1 -uniqueOnly 1 -baq 1 -C 50 -minMapQ 20 -minQ 20 -minInd $((N_IND*1/4)) -setMaxDepth $((N_IND*10)) -doCounts 1 -doGlf 2 -GL 2 -doMajorMinor 4 -doMaf 1 -doPost 2 -doGeno 2 -dumpCounts 2 -postCutoff 0.95 -doHaploCall 1 -out /home/projects/dp_00007/people/hmon/Flat_oysters/02_angsdPopGen/heterozygosity/GEO_FlatOysters--AllSamples_0.25_SITES
```
