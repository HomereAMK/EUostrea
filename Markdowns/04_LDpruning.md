Linkage disequilibrium pruning per chromosome
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
Rscript --vanilla --slave ./BSG_Combined--LD_Server.R
```

```
#variables
REF=/home/projects/dp_00007/people/hmon/AngsdPopStruct/01_infofiles/fileOegenome10scaffoldC3G.fasta
BAMLIST=/home/projects/dp_00007/people/hmon/EUostrea/01_infofiles/bamlist_EUostrea.txt
BASEDIR=/home/projects/dp_00007/people/hmon/EUostrea
OUTPUTFOLDER=/home/projects/dp_00007/people/hmon/EUostrea/03_datasets/LDpruning
N_IND=`cat /home/projects/dp_00007/people/hmon/EUostrea/01_infofiles/bamlist_EUostrea.txt | wc -l`
```


##### Gets .pos file:
```
zcat $BASEDIR/VariantCalling/angsd0.937_htslib1.16_minMapQ20minQ20_minInd145.25_setMinDepthInd5_setMinDepth600setMaxDepth1200.mafs.gz | tail -n +2 | cut -f 1,2 > $OUTPUTFOLDER/angsd0.937_htslib1.16_minMapQ20minQ20_minInd145.25_setMinDepthInd5_setMinDepth600setMaxDepth1200.LD.pos
```

##### Run ngsLD (needs to give the number of sites)
```
/home/projects/dp_00007/apps/ngsLD/ngsLD --n_threads 40 --geno $BASEDIR/VariantCalling/angsd0.937_htslib1.16_minMapQ20minQ20_minInd145.25_setMinDepthInd5_setMinDepth600setMaxDepth1200.beagle.gz --probs --n_ind 582 --n_sites xxx --pos  $OUTPUTFOLDER/angsd0.937_htslib1.16_minMapQ20minQ20_minInd145.25_setMinDepthInd5_setMinDepth600setMaxDepth1200.LD.pos --max_kb_dist 100 | pigz -p 40 >  $OUTPUTFOLDER/angsd0.937_htslib1.16_minMapQ20minQ20_minInd145.25_setMinDepthInd5_setMinDepth600setMaxDepth1200.LD.gz
```

##### Set variables-- the chromosomes
```
CHRs=("scaffold1" "scaffold2" "scaffold3" "scaffold4" "scaffold5" "scaffold6" "scaffold7" "scaffold8" "scaffold9" "scaffold10")
```

##### Splits ngsLD file per chromosome:
```
for query in ${CHRs[*]}
    do
    zcat $OUTPUTFOLDER/angsd0.937_htslib1.16_minMapQ20minQ20_minInd145.25_setMinDepthInd5_setMinDepth600setMaxDepth1200.LD.gz | grep "${query}" | pigz -p 40 > $OUTPUTFOLDER/angsd0.937_htslib1.16_minMapQ20minQ20_minInd145.25_setMinDepthInd5_setMinDepth600setMaxDepth1200.LD.${query}.gz
done
```

##### Get the LD files list
```
find $OUTPUTFOLDER/angsd0.937_htslib1.16_minMapQ20minQ20_minInd145.25_setMinDepthInd5_setMinDepth600setMaxDepth1200.LD.*.gz > $OUTPUTFOLDER/angsd0.937_htslib1.16_minMapQ20minQ20_minInd145.25_setMinDepthInd5_setMinDepth600setMaxDepth1200.LD.PerCHR.list
```

##### LD pruning:
```
parallel --will-cite --dryrun python $BASEDIR/00_scripts/prune_ngsLD.py --input {} --max_dist 100000 --min_weight 0.2 --output  $OUTPUTFOLDER/{/.}.pruneIN --print_excl  $OUTPUTFOLDER/{/.}.pruneOUT :::: $OUTPUTFOLDER/angsd0.937_htslib1.16_minMapQ20minQ20_minInd145.25_setMinDepthInd5_setMinDepth600setMaxDepth1200.LD.PerCHR.list > $OUTPUTFOLDER/RunsLDpruning.txt
```

```
for num in `seq 0 10`
do
num1=`expr $num + 1`
# echo $num1
head -n $num1 $OUTPUTFOLDER/RunsLDpruning.txt | tail -n1 > $OUTPUTFOLDER/SplitLD/RunsLDpruningUp_${num}.sh
done
```

```
for num in `seq 0 10`
do
chmod +x $OUTPUTFOLDER/SplitLD/RunsLDpruningUp_${num}.sh
done
```

##### Merges all pruneIN files:
```
cat /home/projects/dp_00007/people/geopac/Analyses/Turbot/Turbot_LD/BSG_Turbot--AllSamples_Ind30_SNPs.LD.CP*.1.pruneIN > /home/projects/dp_00007/people/geopac/Analyses/Turbot/Turbot_LD/BSG_Turbot--AllSamples_Ind30_SNPs.LD.AllCHRs.pruneIN
```


##### Filters .beagle file based on .pruneIN file:
````
zcat /home/projects/dp_00007/people/geopac/Analyses/Turbot/Turbot_ANGSD/BSG_Turbot--AllSamples_Ind30_SNPs.beagle.gz | awk 'NR==FNR{x[$1]++} NR!=FNR && (FNR == 1 || x[$1]){print}' /home/projects/dp_00007/people/geopac/Analyses/Turbot/Turbot_LD/BSG_Turbot--AllSamples_Ind30_SNPs.LD.AllCHRs.pruneIN - | pigz -p 40 > /home/projects/dp_00007/people/geopac/Analyses/Turbot/Turbot_LD/BSG_Turbot--AllSamples_Ind30_SNPs.Pruned.beagle.gz
````

##### Gets .pos file:
##### Runs ngsLD:
