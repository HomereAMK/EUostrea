Linkage disequilibrium pruning per chromosome TRIAL
================


##### Loads modules:

```
module load tools computerome_utils/2.0
module load perl/5.30.2
module load pigz/2.3.4
module load pandas-profiling/3.0.0
module load graph-tool/2.43
module load gcc/11.1.0
module load intel/perflibs/64/2020_update2
module load R/4.2.0
module load parallel/20220422
#Rscript --vanilla --slave ./BSG_Combined--LD_Server.R
```

```
#variables
REF=/home/projects/dp_00007/people/hmon/AngsdPopStruct/01_infofiles/fileOegenome10scaffoldC3G.fasta
BAMLIST=/home/projects/dp_00007/people/hmon/EUostrea/01_infofiles/bamlist_EUostrea.txt
BASEDIR=/home/projects/dp_00007/people/hmon/EUostrea
OUTPUTFOLDER=/home/projects/dp_00007/people/hmon/EUostrea/03_datasets/LDpruning
N_IND=`cat /home/projects/dp_00007/people/hmon/EUostrea/01_infofiles/bamlist_EUostrea.txt | wc -l`
OUTPUTFOLDER2=/home/projects/dp_00007/people/hmon/EUostrea/03_datasets/Trial
TRIAL=/home/projects/dp_00007/people/hmon/EUostrea/03_datasets/VariantCalling/Dec22_A940_minMapQ20minQ20_NOMININD_setMinDepthInd1_setMinDepthInd1_setMinDepth600setMaxDepth1200
```


##### Gets .pos file:
```
zcat $TRIAL.mafs.gz | tail -n +2 | cut -f 1,2 > $TRIAL.LD.pos
```
🤝
##### Run ngsLD (needs to give the number of sites)
```
/home/projects/dp_00007/apps/ngsLD/ngsLD --n_threads 40 --geno $TRIAL.beagle.gz --probs --n_ind 581 --n_sites 5684643 --pos $TRIAL.LD.pos --max_kb_dist 100 | pigz -p 40 >  $TRIAL.max_kb_dist100.LD.gz
```
🤝
##### Set variables-- the chromosomes
```
CHRs=("scaffold1" "scaffold2" "scaffold3" "scaffold4" "scaffold5" "scaffold6" "scaffold7" "scaffold8" "scaffold9" "scaffold10")
```
🤝w
##### Splits ngsLD file per chromosome:
```
for query in ${CHRs[*]}
    do
    zcat $TRIAL.max_kb_dist100.LD.gz | grep "${query}" | pigz -p 40 > $TRIAL.max_kb_dist100.LD."${query}".gz
done
```
🤝
##### Get the LD files list
```
find $TRIAL.max_kb_dist100.LD.*.gz > $TRIAL.max_kb_dist100.LD.PerCHR.list
```
##### LD pruning:
```
parallel --will-cite --dryrun python $BASEDIR/00_scripts/prune_ngsLD.py --input {} --max_dist 100000 --min_weight 0.2 --output  $OUTPUTFOLDER/{/.}.pruneIN --print_excl  $OUTPUTFOLDER/{/.}.pruneOUT :::: $TRIAL.max_kb_dist100.LD.PerCHR.list > $TRIAL.max_kb_dist100.-min_weight0.2.RunsLDpruning.19jan22.txt
```

```
parallel --will-cite --dryrun python $BASEDIR/00_scripts/prune_ngsLD.py --input {} --max_dist 100000 --min_weight 0.2 --output  $OUTPUTFOLDER/{/.}.19jan23.pruneIN --print_excl  $OUTPUTFOLDER/{/.}.19jan23.pruneOUT :::: $TRIAL.max_kb_dist100.LD.PerCHR.list > $TRIAL.max_kb_dist100kb.-min_weight0.2.19jan23.RunsLDpruning.txt
```

🤝
```
for num in `seq 0 10`
do
    num1=`expr $num + 1`
    # echo $num1
    head -n $num1 $TRIAL.max_kb_dist100.-min_weight0.2.RunsLDpruning.txt | tail -n1 > $TRIAL.max_kb_dist100.-min_weight0.2.SplitLD_RunsLDpruningUp_${num}.sh
done
```
```
for num in `seq 0 10`
do
    num1=`expr $num + 1`
    # echo $num1
    head -n $num1 $TRIAL.max_kb_dist100kb.-min_weight0.2.19jan23.RunsLDpruning.txt | tail -n1 > $TRIAL.max_kb_dist100kb.-min_weight0.2.19jan23.SplitLD_RunsLDpruningUp_${num}.sh
done
```

🤝
```
for num in `seq 0 10`
do
    chmod +x $TRIAL.max_kb_dist100kb.-min_weight0.5.SplitLD_RunsLDpruningUp_${num}.sh
done
```

```
for num in `seq 0 10`
do
    chmod +x $TRIAL.max_kb_dist100kb.-min_weight0.2.19jan23.SplitLD_RunsLDpruningUp_${num}.sh
done
```

🤝
##### Run the $TRIAL.SplitLD_RunsLDpruningUp_${num}.sh jobs:
```
for job in #$TRIAL.max_kb_dist100kb.-min_weight0.5.SplitLD_RunsLDpruningUp_0.sh🤝 \
#$TRIAL.max_kb_dist100kb.-min_weight0.5.SplitLD_RunsLDpruningUp_1.sh🤝 \
$TRIAL.max_kb_dist100kb.-min_weight0.5.SplitLD_RunsLDpruningUp_2.sh \
$TRIAL.max_kb_dist100kb.-min_weight0.5.SplitLD_RunsLDpruningUp_3.sh \
$TRIAL.max_kb_dist100kb.-min_weight0.5.SplitLD_RunsLDpruningUp_4.sh \
#$TRIAL.max_kb_dist100kb.-min_weight0.5.SplitLD_RunsLDpruningUp_5.sh \
$TRIAL.max_kb_dist100kb.-min_weight0.5.SplitLD_RunsLDpruningUp_6.sh \
$TRIAL.max_kb_dist100kb.-min_weight0.5.SplitLD_RunsLDpruningUp_7.sh \
$TRIAL.max_kb_dist100kb.-min_weight0.5.SplitLD_RunsLDpruningUp_8.sh \
$TRIAL.max_kb_dist100kb.-min_weight0.5.SplitLD_RunsLDpruningUp_9.sh \
$TRIAL.max_kb_dist100kb.-min_weight0.5.SplitLD_RunsLDpruningUp_10.sh;
do
    $job
done
```
```
#will only launch the next job if the precedent is finished
job1.sh && job2.sh && job3.sh
```
python prune_ngsLD.py --input /home/projects/dp_00007/people/hmon/EUostrea/03_datasets/VariantCalling/Dec22_A940_minMapQ20minQ20_NOMININD_setMinDepthInd1_setMinDepthInd1_setMinDepth600setMaxDepth1200.max_kb_dist100.LD.scaffold1.gz --max_dist 100000 --min_weight 0.2 --output /home/projects/dp_00007/people/hmon/EUostrea/03_datasets/LDpruning/Dec22_A940_minMapQ20minQ20_NOMININD_setMinDepthInd1_setMinDepthInd1_setMinDepth600setMaxDepth1200.max_kb_dist100.LD.scaffold1.pruneIN --print_excl /home/projects/dp_00007/people/hmon/EUostrea/03_datasets/LDpruning/Dec22_A940_minMapQ20minQ20_NOMININD_setMinDepthInd1_setMinDepthInd1_setMinDepth600setMaxDepth1200.max_kb_dist100.LD.scaffold1.pruneOUT

```
#!/bin/bash
#PBS -d /home/projects/dp_00007/people/hmon/EUostrea
#PBS -W group_list=dp_00007 -A dp_00007
#PBS -N LDprune_minweight0.2.19jan23
#PBS -e 98_log_files/LDprune_minweight0.2.19jan23.err
#PBS -o 98_log_files/LDprune_minweight0.2.19jan23.out
#PBS -l nodes=1:ppn=40:thinnode
#PBS -l walltime=800:00:00
#PBS -l mem=180gb
#PBS -m n
#PBS -r n

module load tools computerome_utils/2.0
module load perl/5.30.2
module load pigz/2.3.4
module load pandas-profiling/3.0.0
module load graph-tool/2.43
module load gcc/11.1.0
module load intel/perflibs/64/2020_update2
module load R/4.2.0
module load parallel/20220422

TRIAL=/home/projects/dp_00007/people/hmon/EUostrea/03_datasets/VariantCalling/Dec22_A940_minMapQ20minQ20_NOMININD_setMinDepthInd1_setMinDepthInd1_setMinDepth600setMaxDepth1200

cd /home/projects/dp_00007/people/hmon/EUostrea/03_datasets/VariantCalling
for num in `seq 0 10`
do
    $TRIAL.max_kb_dist100kb.-min_weight0.2.19jan23.SplitLD_RunsLDpruningUp_${num}.sh
done

```
🤝
##### Merges all pruneIN files:
```
cat /home/projects/dp_00007/people/hmon/EUostrea/03_datasets/LDpruning/Dec22*.pruneIN > /home/projects/dp_00007/people/hmon/EUostrea/03_datasets/LDpruning/Dec22.LD.AllCHRs.min_weight0.5.pruneIN
```

| --- | # of Scaffolds | # PruneIN  | # PruneOUT | # --min_weight |
| :---: | :---: | :---: | :---: |:---: |
| CHRs | 01 | xx | 100% | xx |
| CHRs | 02 | xx | 100% | xx |
| CHRs | 03 | xx | 100% | xx |
| CHRs | 04 | xx | 100% | xx |
| CHRs | 05 | xx | 100% | xx |
| CHRs | 06 | xx | 100% | xx |
| CHRs | 07 | xx | 100% | xx |
| CHRs | 08 | xx | 100% | xx |
| CHRs | 09 | 97790 | 255023 | 0.5 |
| CHRs | 10 | xx | 100% | xx |


From 5684643  sites after LD calculation min weight 0.5: 1404180 sites

scp .pruneIN file + mafs.gz file and run the Rscript Desktop/Scripts/Flat_oysters/04_local_R/00_scripts/NIC_FormattingSnpListforAngsd_11jan22.R with it.


🤝


```R
rm(list=ls(all.names = TRUE))

library(tidyverse)
basedir="/home/projects/dp_00007/people/hmon/EUostrea"  

chr<- c("scaffold1:", "scaffold2:", "scaffold3:", "scaffold4:", "scaffold5:",
        "scaffold6:", "scaffold7:", "scaffold8:", "scaffold9:", "scaffold10:")

pruned_position <- read_lines(paste0(basedir, "/03_datasets/LDpruning/Dec22.LD.AllCHRs.min_weight0.5.pruneIN")) %>% 
  str_remove("scaffold10:") %>% str_remove("scaffold1:") %>% str_remove("scaffold2:") %>% 
  str_remove("scaffold3:") %>% str_remove("scaffold4:") %>% str_remove("scaffold5:") %>%
  str_remove("scaffold6:") %>% str_remove("scaffold7:") %>% str_remove("scaffold8:") %>%
  str_remove("scaffold9:") %>% 
  as.integer()
pruned_position


pruned_snp_list <- read_tsv(paste0(basedir, "/03_datasets/SetAngsdFilters/Jan23_A940_minMapQ20minQ20_NOMININD_setMinDepthInd1_setMinDepthInd1_setMinDepth600setMaxDepth1200.mafs.gz"))%>%
  dplyr::select(1:4) %>%
  filter(position %in% pruned_position)

write_tsv(pruned_snp_list, paste0(basedir, "/03_datasets/LDpruning/LDprunedlist_rightmafs_AllCHRs.min_weight0.5_23jan23"), col_names = F)
```
From 5684643 sites after LD calculation min weight 0.5: 1404180 sites

##### Re-run angsd with the produced Pruned SNPs list with no linked SNPs:



