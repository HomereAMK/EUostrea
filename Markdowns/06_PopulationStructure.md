Population Structure Analysis
================

- [Population Structure Analysis](#population-structure-analysis)
  - [Modules](#modules)
  - [Re-run angsd LD Pruned SNPs list minweight0.5 Global dataset](#re-run-angsd-ld-pruned-snps-list-minweight05-global-dataset)
  - [Missing Data on the LD pruned SNPs](#missing-data-on-the-ld-pruned-snps)
        - [Gets Real Coverage (_Genotype Likelihoods_):](#gets-real-coverage-genotype-likelihoods)
        - [Gets Missing Data (_Genotype Likelihoods_) 11 560 052 SNPs:](#gets-missing-data-genotype-likelihoods-11-560-052-snps)
  - [Re-run angsd LD Pruned SNPs list minweight0.5 Scandinavia](#re-run-angsd-ld-pruned-snps-list-minweight05-scandinavia)
  - [NGSadmix](#ngsadmix)
    - [qsub job because it takes way too much time](#qsub-job-because-it-takes-way-too-much-time)
    - [NGSadmix global dataset](#ngsadmix-global-dataset)
    - [NGSadmix Scandinavia](#ngsadmix-scandinavia)
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
## Re-run angsd LD Pruned SNPs list minweight0.5 Global dataset 
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
> The Global PCA for the MS is coming from this input .covmat...
🤝

## Missing Data on the LD pruned SNPs
##### Gets Real Coverage (_Genotype Likelihoods_):
```bash
angsd -b $BAMLIST -ref $REF -out $OUTPUTFOLDER/30jan23_prunedLDminweight0.5_PopStruct \
-doMajorMinor 3 -doCounts 1 -dumpCounts 2 \
-minQ 20 -minMapQ 20 \
-GL 1 -doGlf 2 \
-P $THREADS \
$EXTRA_ARG \
-sites $SNP_LIST \
-rf $LG_LIST

```
🤝
```bash
zcat $OUTPUTFOLDER/30jan23_prunedLDminweight0.5_PopStruct.counts.gz | tail -n +2 | gawk ' {for (i=1;i<=NF;i++){a[i]+=$i;++count[i]}} END{ for(i=1;i<=NF;i++){print a[i]/count[i]}}' | paste /home/projects/dp_00007/people/hmon/EUostrea/01_infofiles/bamlist_EUostrea.labels - > $OUTPUTFOLDER/30jan23_prunedLDminweight0.5_PopStruct.GL-RealCoverage.txt
```
##### Gets Missing Data (_Genotype Likelihoods_) 11 560 052 SNPs:
```bash
N_SITES=`zcat /home/projects/dp_00007/people/hmon/EUostrea/03_datasets/PopulationStructure/30jan23_prunedLDminweight0.5_PopStruct.beagle.gz | tail -n +2 | wc -l`

zcat /home/projects/dp_00007/people/hmon/EUostrea/03_datasets/PopulationStructure/30jan23_prunedLDminweight0.5_PopStruct.beagle.gz | tail -n +2 | perl /home/projects/dp_00007/apps/Scripts/call_geno.pl --skip 3 | cut -f 4- | awk '{ for(i=1;i<=NF; i++){ if($i==-1)x[i]++} } END{ for(i=1;i<=NF; i++) print i"\t"x[i] }' | paste  /home/projects/dp_00007/people/hmon/EUostrea/01_infofiles/bamlist_EUostrea.2.labels - | awk -v N_SITESawk="$N_SITES" '{print $1"\t"$3"\t"$3*100/N_SITESawk}' | awk '{split($1,b,"-"); print b[1]"\t"$0}' > /home/projects/dp_00007/people/hmon/EUostrea/03_datasets/PopulationStructure/30jan23_prunedLDminweight0.5_PopStruct.GL-MissingData.txt
```


## Re-run angsd LD Pruned SNPs list minweight0.5 Scandinavia
```bash
BAMLISTSCAND=/home/projects/dp_00007/people/hmon/EUostrea/01_infofiles/bamlist_EUostrea_Scandinavia.txt

#Get the beagle file
angsd -b $BAMLISTSCAND -ref $REF -out $OUTPUTFOLDER/2mar23_prunedLDminweight0.5snps_SCAND \
-doMajorMinor 3 -doCounts 1 -doIBS 1 -makematrix 1 -doCov 1 \
-minQ 20 -minMapQ 20 \
-GL 1 -doGlf 2 \
-P $THREADS \
$EXTRA_ARG \
-sites $SNP_LIST \
-rf $LG_LIST
```

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

### NGSadmix global dataset
```bash

BEAGLE=/home/projects/dp_00007/people/hmon/EUostrea/03_datasets/PopulationStructure/30jan23_prunedLDminweight0.5_PopStruct.beagle.gz

for k in $(seq 1 10); do
  /home/projects/dp_00007/apps/Scripts/wrapper_ngsAdmix.sh -P 40 -debug 1 -likes $BEAGLE -K $k -minMaf 0 -tol 1e-6 -tolLike50 1e-3 -maxiter 10000 -o $OUTPUTFOLDER/30jan23_prunedLDminweight0.5_NGSadmix.$k
done

for k in $(seq 8 10); do
  /home/projects/dp_00007/apps/Scripts/wrapper_ngsAdmix.sh -P 40 -debug 1 -likes $BEAGLE -K $k -minMaf 0 -tol 1e-6 -tolLike50 1e-3 -maxiter 2000 -o $OUTPUTFOLDER/27feb23_prunedLDminweight0.5_2Kiter_NGSadmix.$k
done

for k in $(seq 4 10); do
  /home/projects/dp_00007/apps/Scripts/wrapper_ngsAdmix.sh -P 40 -debug 1 -likes $BEAGLE -K $k -minMaf 0 -tol 1e-6 -tolLike50 1e-3 -maxiter 10000 -o $OUTPUTFOLDER/30jan23_prunedLDminweight0.5_NGSadmix.$k
done

```

### NGSadmix Scandinavia
```bash
BEAGLESCAND=/home/projects/dp_00007/people/hmon/EUostrea/03_datasets/PopulationStructure/1mar23_prunedLDminweight0.5snps_SCAND.beagle.gz
for k in $(seq 1 10); do
  /home/projects/dp_00007/apps/Scripts/wrapper_ngsAdmix.sh -P 40 -debug 1 -likes $BEAGLE -K $k -minMaf 0 -tol 1e-6 -tolLike50 1e-3 -maxiter 2000 -o $OUTPUTFOLDER/1mar23_prunedLDminweight0.5snps_SCAND_NGSadmix2kiter.$k
done
```
## evalAdmix
```bash
PREFIXLDPRUNED=/home/projects/dp_00007/people/hmon/EUostrea/03_datasets/PopulationStructure/30jan23_prunedLDminweight0.5_NGSadmix
BEAGLE=/home/projects/dp_00007/people/hmon/EUostrea/03_datasets/PopulationStructure/30jan23_prunedLDminweight0.5_PopStruct.beagle.gz

for k in $(seq 1 20); do
evalAdmix -beagle $BEAGLE -fname $PREFIXLDPRUNED.$k.fopt.gz -qname $PREFIXLDPRUNED.$k.qopt -minMaf 0 -o $PREFIXLDPRUNED.evaladmixOut.$k.corres -P 40
```

```bash
#test evaladmix
evalAdmix -beagle $BEAGLE -fname 1mar23_prunedLDminweight0.5snps_SCAND_NGSadmix2kiter.2.fopt.gz -qname 1mar23_prunedLDminweight0.5snps_SCAND_NGSadmix2kiter.2.qopt -minMaf 0 -o 1mar23_prunedLDminweight0.5snps_SCAND_NGSadmix2kiter.2.corres -P 30
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
🤝