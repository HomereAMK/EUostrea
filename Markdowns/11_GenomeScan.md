Genome scan 
================
  - [Modules](#modules)
  - [Run PCangsd LD Pruned SNPs list minweight0.5 and minMaf 0.05](#run-pcangsd-ld-pruned-snps-list-minweight05-and-minmaf-005)
  - [Run PCangsd with the global Variant calling SNP list](#Run-PCangsd-with-the-global-Variant-calling-SNP-list)

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

```

OUTPUTFOLDER=/home/projects/dp_00007/people/hmon/EUostrea/03_datasets/GenomeScan


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

genome_selection <- npyLoad(paste0(basedir, "/1feb23_GlobalVariants_selectionscan__minmaf0.05_pcangsd_e10.selection.npy"))  

genome_selection_df <- as.data.frame(genome_selection)

basedir="/home/projects/dp_00007/people/hmon/EUostrea/03_datasets/GenomeScan"  

write_tsv(genome_selection_df, paste0(basedir, "/1feb23_GlobalVariants_selectionscan_pcangsd_e10.selection.tsv"), col_names = F)
```
