Fst-based
================
- [Fst-based](#fst-based)
  - [Modules](#modules)
  - [Get minor allele frequencies per population with minIND= 25% and Global SNP list](#get-minor-allele-frequencies-per-population-with-minind-25-and-global-snp-list)
  - [Estimate Fst in each pair of populations](#estimate-fst-in-each-pair-of-populations)
  - [Get “effective sample size” per sample per site](#get-effective-sample-size-per-sample-per-site)
  - [Biological cluster Fst Manhattan plot](#biological-cluster-fst-manhattan-plot)



## Modules
```bash
#qsub jobs info
#!/bin/bash
#PBS -d /home/projects/dp_00007/people/hmon/EUostrea
#PBS -W group_list=dp_00007 -A dp_00007
#PBS -N 09_fstbased.1feb23
#PBS -e 09_fstbased.1feb23.err
#PBS -o 09_fstbased.1feb23.out
#PBS -l nodes=1:ppn=13:thinnode
#PBS -l walltime=400:00:00
#PBS -l mem=100gb
#PBS -m n
#PBS -r n

#angsd
module load tools computerome_utils/2.0
module load htslib/1.16
module load bedtools/2.30.0
module load pigz/2.3.4
module load parallel/20210722
module load angsd/0.940

#variables
REF=/home/projects/dp_00007/people/hmon/AngsdPopStruct/01_infofiles/fileOegenome10scaffoldC3G.fasta
BAMLIST=/home/projects/dp_00007/people/hmon/EUostrea/01_infofiles/bamlist_EUostrea.txt
BASEDIR=/home/projects/dp_00007/people/hmon/EUostrea
OUTPUTFOLDER=/home/projects/dp_00007/people/hmon/EUostrea/03_datasets/Fst
THREADS=10
EXTRA_ARG='-remove_bads 1 -only_proper_pairs 1 -C 50'
LG_LIST=/home/projects/dp_00007/people/hmon/EUostrea/01_infofiles/List_scaffold_28jan23.txt

```

## Get minor allele frequencies per population with minIND= 25% and Global SNP list
> Generate a list of individuals per pops
```bash
POP=("MOLU" "ZECE" "CRES" "ORIS" "CORS" "PONT"  "RIAE" "MORL" "USAM" "TOLL" "COLN" "BARR" "TRAL" "CLEW" "RYAN" "GREV" "WADD" "NISS" "LOGS" "VENO" "HALS" "THIS" "KALV" "HYPP" "LANG" "BUNN" "DOLV" "HAUG" "HAFR" "INNE" "VAGS" "AGAB" "OSTR")

#POP=("INNE" "VAGS")
for query in ${POP[*]}
do 
    grep ${query} $BAMLIST > $BASEDIR/01_infofiles/Jan23--EUostrea_${query}-Fst.list
done
```

> Create a SNP list from the Global Variant calling
```bash
gunzip -c $BASEDIR/03_datasets/SetAngsdFilters/Jan23_A940_minMapQ20minQ20_NOMININD_setMinDepthInd1_setMinDepthInd1_setMinDepth600setMaxDepth1200.mafs.gz | cut -f 1,2,3,4 | tail -n +2 > $OUTPUTFOLDER/global_snp_list_setMinDepth600setMaxDepth1200_jan23.txt
angsd sites index $OUTPUTFOLDER/global_snp_list_setMinDepth600setMaxDepth1200_jan23.txt
```

>This script is used to get minor allele frequency estimation from angsd for each population / group
```bash

for i1 in `seq 0 $((${#POP[@]}-2))`
do
    N_IND=`cat /$BASEDIR/01_infofiles/Jan23--EUostrea_${POP[i1]}-Fst.list | wc -l`
    angsd -nThreads 10 \
    -bam $BASEDIR/01_infofiles/Jan23--EUostrea_${POP[i1]}-Fst.list \
    -anc $REF \
    -ref $REF \
    -out /home/projects/dp_00007/data/hmon/angsd_Fst/EUostrea/31Jan23--mindInd0.25_Unfolded_EUostrea_globalList_${POP[i1]} \
    -dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 3 -doCounts 1 -doDepth 1 -dumpCounts 1 \
    -P $THREADS \
    -minInd $((N_IND*1/4)) -minQ 20 -minMapQ 20 \
    -sites $OUTPUTFOLDER/global_snp_list_setMinDepth600setMaxDepth1200_jan23.txt -rf $LG_LIST \
    $EXTRA_ARG
done
```
🤝 
## Estimate Fst in each pair of populations

> Get the sfs step
```bash
#please re-enter the POP variables
POP=("MOLU" "ZECE" "CRES" "ORIS" "CORS" "PONT"  "RIAE" "MORL" "USAM" "TOLL" "COLN" "BARR" "TRAL" "CLEW" "RYAN" "GREV" "WADD" "NISS" "LOGS" "VENO" "HALS" "THIS" "KALV" "HYPP" "LANG" "BUNN" "DOLV" "HAUG" "HAFR" "INNE" "VAGS" "AGAB" "OSTR")

cd /home/projects/dp_00007/data/hmon/angsd_Fst/EUostrea
for i1 in `seq 0 $((${#POP[@]}-2))`
do
    for i2 in `seq $((i1+1)) $((${#POP[@]}-1))`
    do
        pop1="31Jan23--mindInd0.25_Unfolded_EUostrea_globalList_${POP[i1]}"
        pop2="31Jan23--mindInd0.25_Unfolded_EUostrea_globalList_${POP[i2]}"
        N_SITES=`realSFS print $pop1.saf.idx $pop2.saf.idx | wc -l`
         echo -ne "${POP[i1]}\t${POP[i2]}\t$N_SITES\t"
        if [[ $N_SITES == 0 ]]; then
            echo "NA"
        else
            realSFS $pop1.saf.idx $pop2.saf.idx -fold 1 -P 40 > /home/projects/dp_00007/data/hmon/angsd_Fst/EUostrea/6feb23--mindInd0.25_Unfolded_EUostrea_globalList_${POP[i1]}.${POP[i2]}.sfs
            realSFS fst index $pop1.saf.idx $pop2.saf.idx -sfs /home/projects/dp_00007/data/hmon/angsd_Fst/EUostrea/6feb23--mindInd0.25_Unfolded_EUostrea_globalList_${POP[i1]}.${POP[i2]}.sfs -fold 1 -P 40 -fstout /home/projects/dp_00007/data/hmon/angsd_Fst/EUostrea/6feb23--mindInd0.25_Unfolded_EUostrea_globalList_${POP[i1]}.${POP[i2]}
            realSFS fst stats /home/projects/dp_00007/data/hmon/angsd_Fst/EUostrea/6feb23--mindInd0.25_Unfolded_EUostrea_globalList_${POP[i1]}.${POP[i2]}.fst.idx -P 40
            fi
        done
done > /home/projects/dp_00007/people/hmon/EUostrea/03_datasets/Fst/10feb23--mindInd0.25_globalList_Unfolded_Fst.tsv
```

> Sliding window global list 15kb 15kb
```bash

#please re-enter the POP variables
POP=("MOLU" "ZECE" "CRES" "ORIS" "CORS" "PONT"  "RIAE" "MORL" "USAM" "TOLL" "COLN" "BARR" "TRAL" "CLEW" "RYAN" "GREV" "WADD" "NISS" "LOGS" "VENO" "HALS" "THIS" "KALV" "HYPP" "LANG" "BUNN" "DOLV" "HAUG" "HAFR" "INNE" "VAGS" "AGAB" "OSTR")

for i1 in `seq 0 $((${#POP[@]}-2))`
do
    for i2 in `seq $((i1+1)) $((${#POP[@]}-1))`
     do
        pop1="${POP[i1]}"
        pop2="${POP[i2]}"
        realSFS fst stats2 /home/projects/dp_00007/data/hmon/angsd_Fst/EUostrea/6feb23--mindInd0.25_Unfolded_EUostrea_globalList_${POP[i1]}.${POP[i2]}.fst.idx -win 15000 -step 15000 | cut -f 2- | tail -n +2 | awk '{print $1"\t"$1":"$2"\t"$2-15000"\t"$2"\t"$3"\t"$4}' > /home/projects/dp_00007/data/hmon/angsd_Fst/EUostrea/10feb23--mindInd0.25_Unfolded_EUostrea_globalList_${POP[i1]}.${POP[i2]}_15KB_15KB--Fst.tsv 
    done
done
```


## Get “effective sample size” per sample per site
> Get per-individual depth count at reference filtered SNP list
```bash
angsd -b $BAMLIST \
-anc $REF \
-ref $REF \
-out /home/projects/dp_00007/people/hmon/EUostrea/03_datasets/Depthperind/global_snp_list_setMinDepth600setMaxDepth1200_jan23_ind_depth \
-doCounts 1 -dumpCounts 2 \
-P 8 \
-minMapQ 20 -minQ 20 \
-sites $OUTPUTFOLDER/global_snp_list_setMinDepth600setMaxDepth1200_jan23.txt -rf LG_LIST \
$EXTRA_ARG \
```

>Calculate per-position effective sample size
```R
snp_position <- read_tsv("/home/projects/dp_00007/people/hmon/EUostrea/03_datasets/Depthperind/global_snp_list_setMinDepth600setMaxDepth1200_jan23_ind_depth.pos.gz")
annot <- read_tsv("../sample_lists/sample_table_merged_mincov_contamination_filtered.tsv")
ind_depth <- read_tsv("/home/projects/dp_00007/people/hmon/EUostrea/03_datasets/Depthperind/global_snp_list_setMinDepth600setMaxDepth1200_jan23_ind_depth.counts.gz", skip = 1, col_names = annot$V1) %>%
  bind_cols(snp_position, .) %>%
  dplyr::select(-totDepth, -X155) %>%
  tidyfast::dt_pivot_longer(cols = 3:156, names_to = "V1", values_to = "depth") %>%
  mutate(sample_size=1-0.5^depth)
#write_tsv(ind_depth, "../angsd/bam_list_realigned_mincov_contamination_refbias_filtered_mindp184_maxdp404_minind77_minq20_ind_effective_sample_size.tsv.gz")
pop_effective_sample_size <- ind_depth %>%
  left_join(annot %>% dplyr::select(sample_id_corrected, group), by=c("individual"="V1")) %>%
  group_by(chr, pos, group) %>%
  summarize(sample_size=sum(sample_size))
write_tsv(pop_effective_sample_size, "../angsd/bam_list_realigned_mincov_contamination_refbias_filtered_mindp184_maxdp404_minind77_minq20_pop_effective_sample_size.tsv.gz")
```


## Biological cluster Fst Manhattan plot
>Make the list for biological clusters ADRI, MEDI, ATLA, CHAN, SCAN, NORW, OSTR.
```bash
POP=("MOLU" "ZECE" "CRES" "ORIS" "CORS" "PONT" "RIAE" "MORL" "USAM" "TOLL" "COLN" "BARR" "TRAL" "CLEW" "RYAN" "GREV" "WADD" "NISS" "LOGS" "VENO" "HALS" "THIS" "KALV" "HYPP" "LANG" "BUNN" "DOLV" "HAUG" "HAFR" "INNE" "VAGS" "AGAB" "OSTR")

ADRI=("MOLU" "ZECE" "CRES")
MEDI=("ORIS" "CORS")
ATLA=("PONT" "RIAE")
CHAN=("MORL" "TOLL" "COLN" "BARR" "TRAL" "CLEW" "RYAN")
SCAN=("NISS" "LOGS" "VENO" "HALS" "THIS" "KALV" "HYPP" "LANG" "BUNN" "DOLV" "HAUG" "HAFR")
NORW=("INNE" "VAGS" "AGAB")
OSTR=("OSTR")
for query in ${ADRI[*]}
do 
   cat $BASEDIR/01_infofiles/Jan23--EUostrea_${query}-Fst.list >> $BASEDIR/01_infofiles/Mar23--EUostrea--ADRI-Fst.list
done
for query in ${MEDI[*]}
do 
   cat $BASEDIR/01_infofiles/Jan23--EUostrea_${query}-Fst.list >> $BASEDIR/01_infofiles/Mar23--EUostrea--MEDI-Fst.list
done
for query in ${ATLA[*]}
do 
   cat $BASEDIR/01_infofiles/Jan23--EUostrea_${query}-Fst.list >> $BASEDIR/01_infofiles/Mar23--EUostrea--ATLA-Fst.list
done
for query in ${CHAN[*]}
do 
   cat $BASEDIR/01_infofiles/Jan23--EUostrea_${query}-Fst.list >> $BASEDIR/01_infofiles/Mar23--EUostrea--CHAN-Fst.list
done
for query in ${SCAN[*]}
do 
   cat $BASEDIR/01_infofiles/Jan23--EUostrea_${query}-Fst.list >> $BASEDIR/01_infofiles/Mar23--EUostrea--SCAN-Fst.list
done
for query in ${NORW[*]}
do 
   cat $BASEDIR/01_infofiles/Jan23--EUostrea_${query}-Fst.list >> $BASEDIR/01_infofiles/Mar23--EUostrea--NORW-Fst.list
done
for query in ${OSTR[*]}
do 
   cat $BASEDIR/01_infofiles/Jan23--EUostrea_${query}-Fst.list >> $BASEDIR/01_infofiles/Mar23--EUostrea--OSTR-Fst.list
done
```
>This script is used to get minor allele frequency estimation from angsd for each Biological cluster.
```bash
CLUSTER=("ADRI" "MEDI" "ATLA" "CHAN" "SCAN" "NORW" "OSTR")
for query in ${CLUSTER[*]}
do
    N_IND=`cat $BASEDIR/01_infofiles/Mar23--EUostrea--${query}-Fst.list | wc -l`
    angsd -nThreads 10 \
    -bam $BASEDIR/01_infofiles/Mar23--EUostrea--${query}-Fst.list \
    -anc $REF \
    -ref $REF \
    -out /home/projects/dp_00007/data/hmon/angsd_Fst/EUostrea/Mar23--mindInd0.25_Unfolded_EUostrea_globalList_${query} \
    -dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 3 -doCounts 1 -doDepth 1 -dumpCounts 1 \
    -P $THREADS \
    -minInd $((N_IND*1/4)) -minQ 20 -minMapQ 20 \
    -sites $OUTPUTFOLDER/global_snp_list_setMinDepth600setMaxDepth1200_jan23.txt -rf $LG_LIST \
    $EXTRA_ARG
done
```
🤝
> SFS pairwise-biological-cluster step.
```bash
CLUSTER=("ADRI" "MEDI" "ATLA" "CHAN" "SCAN" "NORW" "OSTR")
cd /home/projects/dp_00007/data/hmon/angsd_Fst/EUostrea
for i1 in `seq 0 $((${#CLUSTER[@]}-2))`
do
    for i2 in `seq $((i1+1)) $((${#CLUSTER[@]}-1))`
    do
        pop1="Mar23--mindInd0.25_Unfolded_EUostrea_globalList_${CLUSTER[i1]}"
        pop2="Mar23--mindInd0.25_Unfolded_EUostrea_globalList_${CLUSTER[i2]}"
        N_SITES=`realSFS print $pop1.saf.idx $pop2.saf.idx | wc -l`
        output_file="/home/projects/dp_00007/data/hmon/angsd_Fst/EUostrea/Mar23--mindInd0.25_Unfolded_EUostrea_globalList_${CLUSTER[i1]}.${CLUSTER[i2]}.sfs"
            if [[ $N_SITES == 0 ]]; then
                echo "NA"
            else
                realSFS $pop1.saf.idx $pop2.saf.idx -fold 1 -P 40 > $output_file
                
                realSFS fst index $pop1.saf.idx $pop2.saf.idx -sfs $output_file -fold 1 -P 40 -fstout /home/projects/dp_00007/data/hmon/angsd_Fst/EUostrea/Mar23--mindInd0.25_Unfolded_EUostrea_globalList_${CLUSTER[i1]}.${CLUSTER[i2]}
                realSFS fst stats /home/projects/dp_00007/data/hmon/angsd_Fst/EUostrea/Mar23--mindInd0.25_Unfolded_EUostrea_globalList_${CLUSTER[i1]}.${CLUSTER[i2]}.fst.idx -P 40
        fi
    done
done > /home/projects/dp_00007/people/hmon/EUostrea/03_datasets/Fst/Mar23_BioCluster--Fstvalues_2.tsv
```

>Fst Sliding window global list 15kb 15kb for pairwise biological cluster 
```bash
CLUSTER=("ADRI" "MEDI" "ATLA" "CHAN" "SCAN" "NORW" "OSTR")
cd /home/projects/dp_00007/data/hmon/angsd_Fst/EUostrea

for i1 in `seq 0 $((${#CLUSTER[@]}-2))`
do
    for i2 in `seq $((i1+1)) $((${#CLUSTER[@]}-1))`
     do
        pop1="${CLUSTER[i1]}"
        pop2="${CLUSTER[i2]}"
        realSFS fst stats2 /home/projects/dp_00007/data/hmon/angsd_Fst/EUostrea/Mar23--mindInd0.25_Unfolded_EUostrea_globalList_${CLUSTER[i1]}.${CLUSTER[i2]}.fst.idx -win 15000 -step 15000 | cut -f 2- | tail -n +2 | awk '{print $1"\t"$1":"$2"\t"$2-15000"\t"$2"\t"$3"\t"$4}' > /home/projects/dp_00007/data/hmon/angsd_Fst/EUostrea/Mar23--mindInd0.25_Unfolded_EUostrea_globalList_${CLUSTER[i1]}.${CLUSTER[i2]}_15KB_15KBwin--Fst.tsv 
    done
done


for i1 in `seq 0 $((${#CLUSTER[@]}-2))`
do
    for i2 in `seq $((i1+1)) $((${#CLUSTER[@]}-1))`
     do
        pop1="${CLUSTER[i1]}"
        pop2="${CLUSTER[i2]}"
        realSFS fst stats2 /home/projects/dp_00007/data/hmon/angsd_Fst/EUostrea/Mar23--mindInd0.25_Unfolded_EUostrea_globalList_${CLUSTER[i1]}.${CLUSTER[i2]}.fst.idx -win 5000 -step 5000 | cut -f 2- | tail -n +2 | awk '{print $1"\t"$1":"$2"\t"$2-5000"\t"$2"\t"$3"\t"$4}' > /home/projects/dp_00007/data/hmon/angsd_Fst/EUostrea/Mar23--mindInd0.25_Unfolded_EUostrea_globalList_${CLUSTER[i1]}.${CLUSTER[i2]}_5KB_5KBwin--Fst.tsv 
    done
done

for i1 in `seq 0 $((${#CLUSTER[@]}-2))`
do
    for i2 in `seq $((i1+1)) $((${#CLUSTER[@]}-1))`
     do
        pop1="${CLUSTER[i1]}"
        pop2="${CLUSTER[i2]}"
        realSFS fst stats2 /home/projects/dp_00007/data/hmon/angsd_Fst/EUostrea/Mar23--mindInd0.25_Unfolded_EUostrea_globalList_${CLUSTER[i1]}.${CLUSTER[i2]}.fst.idx -win 1000 -step 1000 | cut -f 2- | tail -n +2 | awk '{print $1"\t"$1":"$2"\t"$2-1000"\t"$2"\t"$3"\t"$4}' > /home/projects/dp_00007/data/hmon/angsd_Fst/EUostrea/Mar23--mindInd0.25_Unfolded_EUostrea_globalList_${CLUSTER[i1]}.${CLUSTER[i2]}_1KB_1KBwin--Fst.tsv 
    done
done
