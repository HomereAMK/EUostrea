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
🤝
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
🤝
##### LD pruning:
```
parallel --will-cite --dryrun python $BASEDIR/00_scripts/prune_ngsLD.py --input {} --max_dist 100000 --min_weight 0.2 --output  $OUTPUTFOLDER/{/.}.pruneIN --print_excl  $OUTPUTFOLDER/{/.}.pruneOUT :::: $TRIAL.max_kb_dist100.LD.PerCHR.list > $TRIAL.max_kb_dist100.-min_weight0.2.RunsLDpruning.txt
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
🤝
```
for num in `seq 0 10`
do
    chmod +x $TRIAL.max_kb_dist100.-min_weight0.2.SplitLD_RunsLDpruningUp_${num}.sh
done
```
🤝
##### Run the $TRIAL.SplitLD_RunsLDpruningUp_${num}.sh jobs:
```
for job in $TRIAL.max_kb_dist100.-min_weight0.2.SplitLD_RunsLDpruningUp_0.sh $TRIAL.max_kb_dist100.-min_weight0.2.SplitLD_RunsLDpruningUp_1.sh $TRIAL.max_kb_dist100.-min_weight0.2.SplitLD_RunsLDpruningUp_2.sh $TRIAL.max_kb_dist100.-min_weight0.2.SplitLD_RunsLDpruningUp_3.sh $TRIAL.max_kb_dist100.-min_weight0.2.SplitLD_RunsLDpruningUp_4.sh $TRIAL.max_kb_dist100.-min_weight0.2.SplitLD_RunsLDpruningUp_5.sh $TRIAL.max_kb_dist100.-min_weight0.2.SplitLD_RunsLDpruningUp_6.sh $TRIAL.max_kb_dist100.-min_weight0.2.SplitLD_RunsLDpruningUp_7.sh $TRIAL.max_kb_dist100.-min_weight0.2.SplitLD_RunsLDpruningUp_8.sh $TRIAL.max_kb_dist100.-min_weight0.2.SplitLD_RunsLDpruningUp_9.sh $TRIAL.max_kb_dist100.-min_weight0.2.SplitLD_RunsLDpruningUp_10.sh;
do
    $job
done
```
```
$TRIAL.max_kb_dist100.-min_weight0.2.SplitLD_RunsLDpruningUp_0.sh && $TRIAL.max_kb_dist100.-min_weight0.2.SplitLD_RunsLDpruningUp_1.sh && $TRIAL.max_kb_dist100.-min_weight0.2.SplitLD_RunsLDpruningUp_2.sh && $TRIAL.max_kb_dist100.-min_weight0.2.SplitLD_RunsLDpruningUp_3.sh && $TRIAL.max_kb_dist100.-min_weight0.2.SplitLD_RunsLDpruningUp_4.sh && $TRIAL.max_kb_dist100.-min_weight0.2.SplitLD_RunsLDpruningUp_5.sh && $TRIAL.max_kb_dist100.-min_weight0.2.SplitLD_RunsLDpruningUp_6.sh && $TRIAL.max_kb_dist100.-min_weight0.2.SplitLD_RunsLDpruningUp_7.sh && $TRIAL.max_kb_dist100.-min_weight0.2.SplitLD_RunsLDpruningUp_8.sh && $TRIAL.max_kb_dist100.-min_weight0.2.SplitLD_RunsLDpruningUp_9.sh && $TRIAL.max_kb_dist100.-min_weight0.2.SplitLD_RunsLDpruningUp_10.sh
```
🤝
##### Merges all pruneIN files:
```
cat /home/projects/dp_00007/people/hmon/EUostrea/03_datasets/LDpruning/Trial*.pruneIN > /home/projects/dp_00007/people/hmon/EUostrea/03_datasets/LDpruning/Trial.LD.AllCHRs.pruneIN
```
From 128911 sites after LD calculation: 22567 sites
scp .pruneIN file + mafs.gz file and run the Rscript Desktop/Scripts/Flat_oysters/04_local_R/00_scripts/NIC_FormattingSnpListforAngsd_11jan22.R with it.
##### Re-run angsd with the produced Pruned SNPs list with no linked SNPs:



##### Gets .pos file:
##### Runs ngsLD:
