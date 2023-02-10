Population Structure Analysis
================

- [Population Structure Analysis](#population-structure-analysis)
  - [Modules](#modules)
  - [Re-run angsd LD Pruned SNPs list minweight0.5](#re-run-angsd-ld-pruned-snps-list-minweight05)
  - [NGSadmix](#ngsadmix)
    - [qsub job because it takes way too much time](#qsub-job-because-it-takes-way-too-much-time)
  - [evalAdmix](#evaladmix)
  - [Run PCangsd LD Pruned SNPs list minweight0.5 and minMaf 0.05](#run-pcangsd-ld-pruned-snps-list-minweight05-and-minmaf-005)
  - [Run PCangsd with the global Variant calling SNP list](#run-pcangsd-with-the-global-variant-calling-snp-list)


## Modules
```bash
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

#evalAdmix
module load evaladmix/20221125

#Load module for R
module load gsl/2.6
module load imagemagick/7.0.10-13
module load gdal/2.2.3
module load geos/3.8.0
module load jags/4.2.0
module load hdf5
module load netcdf
module load boost/1.74.0
module load openssl/1.0.0
module load lapack
module load udunits/2.2.26
module load proj/7.0.0
module load gcc/10.2.0
module load intel/perflibs/64/2020_update2
module load R/4.0.0

# ngsTools +ngsadmix
module load perl/5.30.2         
module load samtools/1.11
module load ngstools/20190624
module load ngsadmix/32


#variables
REF=/home/projects/dp_00007/people/hmon/AngsdPopStruct/01_infofiles/fileOegenome10scaffoldC3G.fasta
BAMLIST=/home/projects/dp_00007/people/hmon/EUostrea/01_infofiles/bamlist_EUostrea.txt
BASEDIR=/home/projects/dp_00007/people/hmon/EUostrea
N_IND=`cat /home/projects/dp_00007/people/hmon/EUostrea/01_infofiles/bamlist_EUostrea.txt | wc -l`
OUTPUTFOLDER=/home/projects/dp_00007/people/hmon/EUostrea/03_datasets/PopulationStructure
THREADS=10
EXTRA_ARG='-remove_bads 1 -only_proper_pairs 1 -C 50'
SNP_LIST=/home/projects/dp_00007/people/hmon/EUostrea/03_datasets/LDpruning/LDprunedlist_rightmafs_AllCHRs.min_weight0.5_23jan23
LG_LIST=/home/projects/dp_00007/people/hmon/EUostrea/01_infofiles/List_scaffold_28jan23.txt
```

#list of chromosomes/LGs/scaffolds for downstream analysis
```
#cut -f1 $SNP_LIST | sort | uniq > $BASEDIR/01_infofiles/List_scaffold_28jan23.txt
LG_LIST=/home/projects/dp_00007/people/hmon/EUostrea/01_infofiles/List_scaffold_28jan23.txt
```
## Re-run angsd LD Pruned SNPs list minweight0.5
```
angsd sites index $SNPLIST
```
```bash
angsd -b $BAMLIST -ref $REF -out $OUTPUTFOLDER/28jan23_prunedLDminweight0.5_PopStruct \
-doMajorMinor 3 -doCounts 1 -doIBS 1 -makematrix 1 -doCov 1 \
-minQ 20 -minMapQ 20 \
-P $THREADS \
$EXTRA_ARG \
-sites $SNP_LIST \
-rf $LG_LIST
```
🤝
```bash
#Get the beagle file
angsd -b $BAMLIST -ref $REF -out $OUTPUTFOLDER/30jan23_prunedLDminweight0.5_PopStruct \
-doMajorMinor 3 -doCounts 1 -doIBS 1 -makematrix 1 -doCov 1 \
-minQ 20 -minMapQ 20 \
-GL 1 -doGlf 2 \
-P $THREADS \
$EXTRA_ARG \
-sites $SNP_LIST \
-rf $LG_LIST
```

🤝



## NGSadmix
### qsub job because it takes way too much time
```bash
#!/bin/bash
#PBS -d /home/projects/dp_00007/people/hmon/EUostrea
#PBS -W group_list=dp_00007 -A dp_00007
#PBS -N NGSAdmix_bignode_7feb23
#PBS -e NGSAdmix_bignode_7feb23.err
#PBS -o NGSAdmix_bignode_7feb23.out
#PBS -l nodes=2:ppn=36:fatnode
#PBS -l walltime=600:00:00
#PBS -l mem=1300gb
#PBS -m n
#PBS -r n
```

```bash

BEAGLE=/home/projects/dp_00007/people/hmon/EUostrea/03_datasets/PopulationStructure/30jan23_prunedLDminweight0.5_PopStruct.beagle.gz

for k in $(seq 1 10); do
  /home/projects/dp_00007/apps/Scripts/wrapper_ngsAdmix.sh -P 40 -debug 1 -likes $BEAGLE -K $k -minMaf 0 -tol 1e-6 -tolLike50 1e-3 -maxiter 10000 -o $OUTPUTFOLDER/30jan23_prunedLDminweight0.5_NGSadmix.$k
done

for k in $(seq 2 4); do
  /home/projects/dp_00007/apps/Scripts/wrapper_ngsAdmix.sh -P 40 -debug 1 -likes $BEAGLE -K $k -minMaf 0 -tol 1e-6 -tolLike50 1e-3 -maxiter 10000 -o $OUTPUTFOLDER/30jan23_prunedLDminweight0.5_NGSadmix.$k
done

for k in $(seq 4 10); do
  /home/projects/dp_00007/apps/Scripts/wrapper_ngsAdmix.sh -P 40 -debug 1 -likes $BEAGLE -K $k -minMaf 0 -tol 1e-6 -tolLike50 1e-3 -maxiter 10000 -o $OUTPUTFOLDER/30jan23_prunedLDminweight0.5_NGSadmix.$k
done

```

## evalAdmix
```bash
LDPRUNED=$OUTPUTFOLDER/30jan23_prunedLDminweight0.5_PopStruct

for k in $(seq 1 20); do
evalAdmix -beagle $LDPRUNED.beagle.gz -fname $LDPRUNED.$k.fopt.gz -qname $LDPRUNED.$k.qopt -o evaladmixOut.$LDPRUNED.$k.corres -P 40
```
## Run PCangsd LD Pruned SNPs list minweight0.5 and minMaf 0.05
```bash
pcangsd -b $OUTPUTFOLDER/30jan23_prunedLDminweight0.5_PopStruct.beagle.gz \
--selection \
--minMaf 0.05 \
--sites_save \
--threads 40 \
-o $OUTPUTFOLDER/30jan23_prunedLDminweight0.5_PopStruct_pcangsd
```
🤝

#on c2
```R
rm(list=ls(all.names = TRUE))

pacman::p_load(vegan, tidyverse, RcppCNPy, pheatmap, extrafont, ggforce, ggrepel, ggstar, np, reticulate, cowplot)


basedir="/home/projects/dp_00007/people/hmon/EUostrea/03_datasets/PopulationStructure"  

genome_selection <- npyLoad(paste0(basedir, "/30jan23_prunedLDminweight0.5_PopStruct_pcangsd.selection.npy"))  

genome_selection_df <- as.data.frame(genome_selection)

write_tsv(genome_selection_df, paste0(basedir, "/30jan23_prunedLDminweight0.5_PopStruct_pcangsd.selection.tsv"), col_names = F)
```
🤝

## Run PCangsd with the global Variant calling SNP list
```bash
GLOBAL_VC_BEAGLE=/home/projects/dp_00007/people/hmon/EUostrea/03_datasets/SetAngsdFilters/Jan23_A940_minMapQ20minQ20_NOMININD_setMinDepthInd1_setMinDepthInd1_setMinDepth600setMaxDepth1200.beagle.gz

pcangsd -b $GLOBAL_VC_BEAGLE \
--selection \
--minMaf 0.05 \
--sites_save \
--threads 40 \
-e 10 \
-o $OUTPUTFOLDER/1feb23_GlobalVariants_selectionscan__minmaf0.05_pcangsd_e10
```

#on c2
```R
rm(list=ls(all.names = TRUE))

pacman::p_load(vegan, tidyverse, RcppCNPy, pheatmap, extrafont, ggforce, ggrepel, ggstar, np, reticulate, cowplot)


basedir="/home/projects/dp_00007/people/hmon/EUostrea/03_datasets/PopulationStructure"  

genome_selection <- npyLoad(paste0(basedir, "/1feb23_GlobalVariants_selectionscan_pcangsd_e10.selection.npy"))  

genome_selection_df <- as.data.frame(genome_selection)

write_tsv(genome_selection_df, paste0(basedir, "/1feb23_GlobalVariants_selectionscan_pcangsd_e10.selection.tsv"), col_names = F)
```
