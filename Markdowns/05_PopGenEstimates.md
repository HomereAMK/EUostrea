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
  - [Heterozygosity](#heterozygosity)
- [Runs \_ANGSD with depth filters to get the mafs file](#runs-_angsd-with-depth-filters-to-get-the-mafs-file)
- [Generates a `.bed` file based on the `.mafs` file:](#generates-a-bed-file-based-on-the-mafs-file)
- [Creates a position file based on this new `.bed`:](#creates-a-position-file-based-on-this-new-bed)
- [Indexs the `.pos` file created above:](#indexs-the-pos-file-created-above)
- [Runs _ANGSD_ under -doSaf::](#runs-angsd-under--dosaf)
- [Calculates fractions:](#calculates-fractions)
- [Calculates heterozygosity:](#calculates-heterozygosity)




## Load module angsd.
```
#angsd
module load tools computerome_utils/2.0
module load htslib/1.16
module load bedtools/2.30.0
module load pigz/2.3.4
module load parallel/20220422
module load angsd/0.940

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
        angsd -nThreads 40 -ref $REF -anc $REF -bam /home/projects/dp_00007/people/hmon/EUostrea/01_infofiles/EUostrea_${query}-Fst.list -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 -minMapQ 20 -minQ 20 -setMinDepthInd 1 -setMinDepth 600 -setMaxDepth 1200 -GL 1 -doSaf 1 -out /home/projects/dp_00007/data/hmon/angsd_PopGen/Jan23--Unfolded_PopGenNIC_${query}
    done
done
```
for query in ${POP[*]}
    do
        REF=/home/projects/dp_00007/people/hmon/AngsdPopStruct/01_infofiles/fileOegenome10scaffoldC3G.fasta
        /home/projects/dp_00007/apps/Scripts/wrapper_angsd.sh -debug 2 -nThreads 40 -ref $REF -anc $REF -bam /home/projects/dp_00007/people/hmon/EUostrea/01_infofiles/EUostrea_${query}-Fst.list -remove_bads 1 -uniqueOnly 1 -baq 1 -C 50 -minMapQ 20 -minQ 20 -minInd $((N_IND*2/3)) -GL 1 -doSaf 1 -out /home/projects/dp_00007/data/hmon/angsd_PopGen/Jan23--Unfolded_PopGenGEO_${query}
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
# Runs _ANGSD with depth filters to get the mafs file
```
N_IND=`cat /home/projects/dp_00007/people/hmon/EUostrea/01_infofiles/bamlist_EUostrea.txt | wc -l`
REF=/home/projects/dp_00007/people/hmon/AngsdPopStruct/01_infofiles/fileOegenome10scaffoldC3G.fasta
BAMLIST=/home/projects/dp_00007/people/hmon/EUostrea/01_infofiles/bamlist_EUostrea.txt

angsd -nThreads 40 -ref $REF -bam $BAMLIST -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
-doCounts 1 -dumpCounts 2 \
-minMapQ 20 -minQ 20 -setMinDepthInd 1 -setMinDepth 600 -setMaxDepth 1200 \
-GL 1 -doMajorMinor 1 -doMaf 1 \
-out /home/projects/dp_00007/people/hmon/EUostrea/03_datasets/Het/Jan23_NICdepthfilters_het_angsd


/home/projects/dp_00007/apps/Scripts/wrapper_angsd.sh -debug 2 -nThreads 40 -ref $REF -bam $BAMLIST -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
-doCounts 1 -dumpCounts 2 \
-minMapQ 20 -minQ 20 -setMinDepthInd 1 -setMinDepth 600 -setMaxDepth 1200 \
-GL 1 -doMajorMinor 1 -doMaf 1 \
-out /home/projects/dp_00007/people/hmon/EUostrea/03_datasets/Het/Jan23_wrap_NICdepthfilters_het_angsd
```
🤝
# Generates a `.bed` file based on the `.mafs` file:
```
zcat /home/projects/dp_00007/people/hmon/EUostrea/03_datasets/Het/Jan23_wrap_NICdepthfilters_het_angsd.mafs.gz | cut -f1,2 | tail -n +2 | awk '{print $1"\t"$2-1"\t"$2}' | bedtools merge -i - > /home/projects/dp_00007/people/hmon/EUostrea/03_datasets/Het/Jan23_wrap_NICdepthfilters_het_angsd.bed
```
🤝
# Creates a position file based on this new `.bed`:
```
awk '{print $1"\t"($2+1)"\t"$3}'  /home/projects/dp_00007/people/hmon/EUostrea/03_datasets/Het/Jan23_wrap_NICdepthfilters_het_angsd.bed > /home/projects/dp_00007/people/hmon/EUostrea/03_datasets/Het/Jan23_wrap_NICdepthfilters_het_angsd_bed.pos
```
🤝
# Indexs the `.pos` file created above:
```
angsd sites index /home/projects/dp_00007/people/hmon/EUostrea/03_datasets/Het/Jan23_wrap_NICdepthfilters_het_angsd_bed.pos
```
🤝
# Runs _ANGSD_ under -doSaf::
```
parallel --plus angsd -i {} -ref $REF -anc $REF -sites /home/projects/dp_00007/people/hmon/EUostrea/03_datasets/Het/Jan23_wrap_NICdepthfilters_het_angsd_bed.pos -GL 1 -doSaf 1 -remove_bads 1 -uniqueOnly 1 -baq 1 -C 50 -minMapQ 20 -minQ 20 -out /home/projects/dp_00007/data/hmon/angsd_Het/{/...} :::: /home/projects/dp_00007/people/hmon/EUostrea/01_infofiles/bamlist_EUostrea.txt
```
🤝
# Calculates fractions:
```
parallel --tmpdir /home/projects/dp_00007/data/hmon/angsd_Het/ -j 1 --plus "realSFS -fold 1 -P 10 {} > /home/projects/dp_00007/data/hmon/angsd_Het/{/..}.het" ::: /home/projects/dp_00007/data/hmon/angsd_Het/*.saf.idx

#alternative 

# Run 40 instances of realSFS in parallel using specified tmpdir and plus options 
parallel --tmpdir /home/projects/dp_00007/data/hmon/angsd_Het/ --plus "realSFS -fold 1 -P 10 {} > /home/projects/dp_00007/data/hmon/angsd_Het/{/..}.het" -j 1 --memfree 160000 ::: /home/projects/dp_00007/data/hmon/angsd_Het/*.saf.idx
```
🤝
# Calculates heterozygosity:
```
fgrep '.' *.23jan23.het | tr ":" " " | awk '{print $1"\t"$3/($2+$3)*100}' | sed -r 's/.het//g' | awk '{split($0,a,"_"); print $1"\t"a[1]"\t"$2"\t"$3'} | awk '{split($2,a,"-"); print $1"\t"a[2]"\t"$3'} > /home/projects/dp_00007/people/hmon/EUostrea/03_datasets/Het/GEO_EUostrea--AllSamples_setMinDepth600_setMaxDepth1200--HET.23jan23.txt
```


