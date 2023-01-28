Population Structure Analysis
================

- [Population Structure Analysis](#population-structure-analysis)
  - [Modules](#modules)
  - [Re-run angsd LD Pruned SNPs list minweight0.5](#re-run-angsd-ld-pruned-snps-list-minweight05)
  - [Run PCangsd LD Pruned SNPs list minweight0.5 and minMaf 0.05](#run-pcangsd-ld-pruned-snps-list-minweight05-and-minmaf-005)


## Modules
```
#angsd
module load tools computerome_utils/2.0
module load htslib/1.16
module load bedtools/2.30.0
module load pigz/2.3.4
module load parallel/20210722
module load angsd/0.940

#pcangsd
module load tools computerome_utils/2.0
module load pcangsd/20220330 
```
#list of chromosomes/LGs/scaffolds for downstream analysis
```
cut -f1 $SNP_LIST | sort | uniq > $BASEDIR/01_infofiles/List_scaffold_28jan23.txt
SNP_LIST=/home/projects/dp_00007/people/hmon/EUostrea/03_datasets/LDpruning/LDprunedlist_rightmafs_AllCHRs.min_weight0.5_23jan23
LG_LIST=/home/projects/dp_00007/people/hmon/EUostrea/01_infofiles/List_scaffold_28jan23.txt
```
## Re-run angsd LD Pruned SNPs list minweight0.5
```
REF=/home/projects/dp_00007/people/hmon/AngsdPopStruct/01_infofiles/fileOegenome10scaffoldC3G.fasta
BAMLIST=/home/projects/dp_00007/people/hmon/EUostrea/01_infofiles/bamlist_EUostrea.txt
BASEDIR=/home/projects/dp_00007/people/hmon/EUostrea
N_IND=`cat /home/projects/dp_00007/people/hmon/EUostrea/01_infofiles/bamlist_EUostrea.txt | wc -l`
OUTPUTFOLDER=/home/projects/dp_00007/people/hmon/EUostrea/03_datasets/PopulationStructure
THREADS=8
EXTRA_ARG='-remove_bads 1 -only_proper_pairs 1 -C 50'
SNP_LIST=/home/projects/dp_00007/people/hmon/EUostrea/03_datasets/LDpruning/LDprunedlist_rightmafs_AllCHRs.min_weight0.5_23jan23
```
```
sites index $SNPLIST
```
```
angsd -b $BAMLIST -ref $REF -out $OUTPUTFOLDER/28jan23_prunedLDminweight0.5_PopStruct \
-doMajorMinor 3 -doCounts 1 -doIBS 1 -makematrix 1 -doCov 1 \
-minQ 20 -minMapQ 20 \
-P $THREADS \
$EXTRA_ARG \
-sites $SNP_LIST \
-rf $LG_LIST \
```
## Run PCangsd LD Pruned SNPs list minweight0.5 and minMaf 0.05
```
pcangsd -b $OUTPUTDIR/28jan23_prunedLDminweight0.5_PopStruct.beagle.gz \
--selection \
--minMaf 0.05 \
--sites_save \
--threads 40 \
-o $OUTPUTDIR/28jan23_PCangsd_prunedLDminweight0.5_PopStruct
```