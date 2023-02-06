Genome scan 
================
- [Genome scan](#genome-scan)
  - [Modules](#modules)
  - [Run PCangsd LD Pruned SNPs list minweight0.5 and minMaf 0.05](#run-pcangsd-ld-pruned-snps-list-minweight05-and-minmaf-005)
  - [Run PCangsd with the global Variant calling SNP list](#run-pcangsd-with-the-global-variant-calling-snp-list)
  - [PCA in inversion-like regions](#pca-in-inversion-like-regions)
    - [Inversions-like regions](#inversions-like-regions)
    - [Get the beagle file](#get-the-beagle-file)

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

#bcftools
module load bcftools/1.16

```



## Run PCangsd LD Pruned SNPs list minweight0.5 and minMaf 0.05
```bash
OUTPUTFOLDER=/home/projects/dp_00007/people/hmon/EUostrea/03_datasets/GenomeScan

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
OUTPUTFOLDER=/home/projects/dp_00007/people/hmon/EUostrea/03_datasets/GenomeScan

GLOBAL_VC_BEAGLE=/home/projects/dp_00007/people/hmon/EUostrea/03_datasets/SetAngsdFilters/Jan23_A940_minMapQ20minQ20_NOMININD_setMinDepthInd1_setMinDepthInd1_setMinDepth600setMaxDepth1200.beagle.gz
GLOBAL_VC_PATH=/home/projects/dp_00007/people/hmon/EUostrea/03_datasets/SetAngsdFilters/Jan23_A940_minMapQ20minQ20_NOMININD_setMinDepthInd1_setMinDepthInd1_setMinDepth600setMaxDepth1200

## Create a SNP list to use in downstream analyses
gunzip -c $GLOBAL_VC_PATH.mafs.gz | cut -f 1,2,3,4 | tail -n +2 > $OUTPUTFOLDER/1feb23_GlobalVariants_selectionscan__minmaf0.05_List.txt

angsd sites index $OUTPUTFOLDER/1feb23_GlobalVariants_selectionscan__minmaf0.05_List.txt

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

genome_selection <- npyLoad(paste0(basedir, "/1feb23_GlobalVariants_selectionscan__minmaf0.05_pcangsd_e10.selection.npy"))  

genome_selection_df <- as.data.frame(genome_selection)

basedir="/home/projects/dp_00007/people/hmon/EUostrea/03_datasets/GenomeScan"  

write_tsv(genome_selection_df, paste0(basedir, "/1feb23_GlobalVariants_selectionscan_pcangsd_e10.selection.tsv"), col_names = F)
```


## PCA in inversion-like regions
> Potential inversion-like regions are visually assessed on the PCAngsd genome scan with -selection flag.
> 3 regions with large number block of SNPs showing signs of selection, with distributions that exceeded what was expected under neutrality along the first principal component: scaffold4 (first 10Mb), scaffold5 (first 20Mb) and scaffold8 (last 30Mb).

### Inversions-like regions
```bash
#Variables
GLOBAL_VC_BEAGLE=/home/projects/dp_00007/people/hmon/EUostrea/03_datasets/SetAngsdFilters/Jan23_A940_minMapQ20minQ20_NOMININD_setMinDepthInd1_setMinDepthInd1_setMinDepth600setMaxDepth1200.beagle.gz
LIST_VC=/home/projects/dp_00007/people/hmon/EUostrea/03_datasets/Fst/global_snp_list_setMinDepth600setMaxDepth1200_jan23.txt
OUTPUTFOLDER=/home/projects/dp_00007/people/hmon/EUostrea/03_datasets/InversionsLike/
```
```bash
# Find the starting and ending position of the first 10 megabases of scaffold 4
lenght_sca4=96428841
start_pos=1
end_pos=$((start_pos + 10**7))

# Extract the region of interest from the global variant calling snp list with awk
cat $LIST_VC | awk -v start=$start_pos -v end=$end_pos '$1 == "scaffold4" && $2 >= start && $2 <= end {print $1,$2,$3,$4}' > $OUTPUTFOLDER/scaffold4/scaffold4_InvReg_awk.txt
```

```bash
# Find the starting and ending position of the first 20 megabases of scaffold 5
lenght_sca5=96132218
start_pos=1
end_pos=$((start_pos + 20**7))

# Extract the region of interest from the global variant calling snp list with awk
cat $LIST_VC | awk -v start=$start_pos -v end=$end_pos '$1 == "scaffold5" && $2 >= start && $2 <= end {print $1,$2,$3,$4}' > $OUTPUTFOLDER/scaffold5/scaffold5_InvReg_awk.txt
```

```bash
# Find the starting position at around 32 megabases of scaffold 8 and the end position is the end of the scaffold 8
lenght_sca8=58502465
start_pos=32000000
end_pos=$lenght_sca8

# Extract the region of interest from the global variant calling snp list with awk
cat $LIST_VC | awk -v start=$start_pos -v end=$end_pos '$1 == "scaffold8" && $2 >= start && $2 <= end {print $1,$2,$3,$4}' > $OUTPUTFOLDER/scaffold8/scaffold8_InvReg_awk.txt
```
### Get the beagle file
```bash
angsd -b $BAMLIST -ref $REF -out $OUTPUTFOLDER/30jan23_prunedLDminweight0.5_PopStruct \
-doMajorMinor 3 -doCounts 1 -doIBS 1 -makematrix 1 -doCov 1 \
-minQ 20 -minMapQ 20 \
-GL 1 -doGlf 2 \
-P $THREADS \
$EXTRA_ARG \
-sites $SNP_LIST \
-rf $LG_LIST
```
