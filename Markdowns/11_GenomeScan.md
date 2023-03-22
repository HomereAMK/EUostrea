Genome scan 
================
- [Genome scan](#genome-scan)
  - [Modules](#modules)
  - [Run PCangsd LD Pruned SNPs list minweight0.5 and minMaf 0.05](#run-pcangsd-ld-pruned-snps-list-minweight05-and-minmaf-005)
  - [Run PCangsd with the global Variant calling SNP list](#run-pcangsd-with-the-global-variant-calling-snp-list)
  - [PCA in inversion-like regions](#pca-in-inversion-like-regions)
    - [Inversions-like regions identification](#inversions-like-regions-identification)
    - [Snps per inversion-like regions](#snps-per-inversion-like-regions)
    - [Get the beagle file](#get-the-beagle-file)
    - [Get the PopGenestimates in the inversion region](#get-the-popgenestimates-in-the-inversion-region)
  - [Outliers on PC9 scaffold2 and 10](#outliers-on-pc9-scaffold2-and-10)
    - [Get the beagle file](#get-the-beagle-file-1)


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

### Inversions-like regions identification
```bash
#Variables
#GLOBAL_VC_BEAGLE=/home/projects/dp_00007/people/hmon/EUostrea/03_datasets/SetAngsdFilters/Jan23_A940_minMapQ20minQ20_NOMININD_setMinDepthInd1_setMinDepthInd1_setMinDepth600setMaxDepth1200.beagle.gz
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
end_pos=20000000
# Extract the region of interest from the global variant calling snp list with awk
cat $LIST_VC | awk -v start=$start_pos -v end=$end_pos '$1 == "scaffold5" && $2 >= start && $2 <= end {print $1,$2,$3,$4}' > $OUTPUTFOLDER/scaffold5/scaffold5_InvReg_awk.txt
```
🤝
```bash
# Find the starting position at around 32 megabases of scaffold 8 and the end position is the end of the scaffold 8
lenght_sca8=58502465
start_pos=32000000
end_pos=$lenght_sca8

# Extract the region of interest from the global variant calling snp list with awk
cat $LIST_VC | awk -v start=$start_pos -v end=$end_pos '$1 == "scaffold8" && $2 >= start && $2 <= end {print $1,$2,$3,$4}' > $OUTPUTFOLDER/scaffold8/scaffold8_InvReg_awk.txt
```
🤝
### Snps per inversion-like regions
| "scaffold4 inv"  | "scaffold5 inv" |"scaffold8 inv" |
| ------------- | ------------- |------------- |
| SNPs  | SNPs  | SNPs  |
| 85239  | 171821  | 191366 |


### Get the beagle file
```bash
#Variables
LG_LIST=/home/projects/dp_00007/people/hmon/EUostrea/01_infofiles/List_scaffold_28jan23.txt
REF=/home/projects/dp_00007/people/hmon/AngsdPopStruct/01_infofiles/fileOegenome10scaffoldC3G.fasta
BAMLIST=/home/projects/dp_00007/people/hmon/EUostrea/01_infofiles/bamlist_EUostrea.txt
THREADS=10
```
🤝
```bash
for query in scaffold8 scaffold4 scaffold5
do
  OUTPUTFOLDER="/home/projects/dp_00007/people/hmon/EUostrea/03_datasets/InversionsLike/${query}"
  SNP_LIST="$OUTPUTFOLDER/${query}_InvReg_awk.txt"
  # Error handling for the 'angsd sites index' command
  if ! angsd sites index "$SNP_LIST"; then
    echo "Error: angsd sites index failed for $query"
    continue
  fi
  # Error handling for the 'angsd' command
  if ! angsd -b "$BAMLIST" -ref "$REF" -out "$OUTPUTFOLDER/6feb23_${query}_InvReg" \
        -doMajorMinor 3 -doCounts 1 -doIBS 1 -makematrix 1 -doCov 1 \
        -minQ 20 -minMapQ 20 \
        -GL 1 -doGlf 2 \
        -P "$THREADS" \
        -sites "$SNP_LIST" \
        -rf /home/projects/dp_00007/people/hmon/Flat_oysters/01_infofiles/${query}.txt ; then
    echo "Error: angsd failed for $query"
  fi
done
```
🤝
### Get the PopGenestimates in the inversion region
>Reg04
scaffold4:14226-9973860
>Reg05
scaffold5:375-19999685
>Reg08
scaffold8:32001789-58488037
```bash
for query in scaffold8 scaffold4 scaffold5
do
  OUTPUTFOLDER="/home/projects/dp_00007/people/hmon/EUostrea/03_datasets/InversionsLike/${query}"

  pcangsd -b "$OUTPUTFOLDER/6feb23_${query}_InvReg.beagle.gz" \
  --selection \
  --minMaf 0.05 \
  --sites_save \
  --threads 40 \
  -e 10 \
  -o "$OUTPUTFOLDER/7feb23_${query}_InvReg_pcangsd"

done
```

#on c2
```R
# Clear the environment
rm(list = ls(all.names = TRUE))

# Load packages
pacman::p_load(vegan, tidyverse, RcppCNPy, pheatmap, extrafont, ggforce, ggrepel, ggstar, np, reticulate, cowplot)

# Define a vector of scaffold names
InvReg <- c("scaffold4", "scaffold8", "scaffold5")

# Loop through the scaffold names
for (query in InvReg) {
  basedir <- paste0("/home/projects/dp_00007/people/hmon/EUostrea/03_datasets/InversionsLike/", query)  

  # Load the npy file
  genome_selection <- RcppCNPy::npyLoad(paste0(basedir, "/7feb23_", query, "_InvReg_pcangsd.selection.npy"))  

  # Convert the npy file to a data frame
  genome_selection_df <- as.data.frame(genome_selection)

  # Write the data frame to a tsv file
  write.table(genome_selection_df, paste0(basedir, "/1feb23_GlobalVariants_InvReg", query, "_pcangsd.selection.tsv"), 
              col.names = FALSE, row.names = FALSE, sep = "\t")
}
```


## Outliers on PC9 scaffold2 and 10
```bash
#Variables
#GLOBAL_VC_BEAGLE=/home/projects/dp_00007/people/hmon/EUostrea/03_datasets/SetAngsdFilters/Jan23_A940_minMapQ20minQ20_NOMININD_setMinDepthInd1_setMinDepthInd1_setMinDepth600setMaxDepth1200.beagle.gz
LIST_VC=/home/projects/dp_00007/people/hmon/EUostrea/03_datasets/Fst/global_snp_list_setMinDepth600setMaxDepth1200_jan23.txt
OUTPUTFOLDER=/home/projects/dp_00007/people/hmon/EUostrea/03_datasets/GenomeScan/Outliers
#outliers PC9 scaffold2
start_pos=93220504
end_pos=93743877
# Extract the region of interest from the global variant calling snp list with awk
cat $LIST_VC | awk -v start=$start_pos -v end=$end_pos '$1 == "scaffold2" && $2 >= start && $2 <= end {print $1,$2,$3,$4}' > $OUTPUTFOLDER/Outliers_PC2_sca2_10_10feb23.txt
#outliers PC9 scaffold10
start_pos=11747526
end_pos=11856952
# Append the region of interest from the global variant calling snp list with awk
cat $LIST_VC | awk -v start=$start_pos -v end=$end_pos '$1 == "scaffold10" && $2 >= start && $2 <= end {print $1,$2,$3,$4}' >> $OUTPUTFOLDER/Outliers_PC2_sca2_10_10feb23.txt

#error on the output names
for file in $OUTPUTFOLDER/Outliers_PC2_sca2_10_10feb23*; do
    mv "$file" "${file/Outliers_PC2/Outliers_PC9}"
done

```
### Get the beagle file
```bash
#Variables
LG_LIST=/home/projects/dp_00007/people/hmon/EUostrea/03_datasets/GenomeScan/Outliers/lg_list_sca2_sca10
REF=/home/projects/dp_00007/people/hmon/AngsdPopStruct/01_infofiles/fileOegenome10scaffoldC3G.fasta
BAMLIST=/home/projects/dp_00007/people/hmon/EUostrea/01_infofiles/bamlist_EUostrea.txt
THREADS=12
OUTLIERS_LIST=/home/projects/dp_00007/people/hmon/EUostrea/03_datasets/GenomeScan/Outliers/Outliers_PC2_sca2_10_10feb23.txt
OUTPUTFOLDER=/home/projects/dp_00007/people/hmon/EUostrea/03_datasets/GenomeScan/Outliers

angsd sites index $OUTLIERS_LIST

angsd -b "$BAMLIST" -ref "$REF" -out "$OUTPUTFOLDER/Outliers_PC2_sca2_10_10feb23" \
-doMajorMinor 3 -doCounts 1 -doIBS 1 -makematrix 1 -doCov 1 \
-minQ 20 -minMapQ 20 \
-GL 1 -doGlf 2 \
-P "$THREADS" \
-sites "$OUTLIERS_LIST" \
-rf $LG_LIST