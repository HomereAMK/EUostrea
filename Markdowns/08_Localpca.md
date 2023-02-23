Local-pca
================

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
```bash

#!/bin/bash
## This script is used to run local_pca based on genotype likelihood data. See https://github.com/petrelharp/local_pca for details. 
## This script will use a separate thread for each LG. So you will need to first run /workdir/genomic-data-analysis/scripts/subset_beagle_by_lg.sh

#variables
BEAGLE=/home/projects/dp_00007/people/hmon/EUostrea/03_datasets/LocalPca/Jan23_A940_minMapQ20minQ20_NOMININD_setMinDepthInd1_setMinDepthInd1_setMinDepth600setMaxDepth1200.beagle.gz
LG_LIST=/home/projects/dp_00007/people/hmon/EUostrea/01_infofiles/List_scaffold_28jan23.txt
SNP=1000 ## Number of SNPs to include in each window
PC=2 ## Number of PCs to keep for each window
N_CORE_MAX=20 # Maximum number of threads to use simulatenously
LOCAL_PCA_1=/home/projects/dp_00007/people/hmon/EUostrea/00_scripts/local_pca_1.sh
LOCAL_PCA_2=/home/projects/dp_00007/people/hmon/EUostrea/00_scripts/local_pca_2.R
PYTHON=python
PCANGSD=pcangsd

## Extract prefix and directory from the beagle path
PREFIX=`echo $BEAGLE | sed 's/\..*//' | sed -e 's#.*/\(\)#\1#'`
BEAGLEDIR=`echo $BEAGLE | sed 's:/[^/]*$::' | awk '$1=$1"/"'`
🤝

COUNT=0
for LG in `cat $LG_LIST`; do
	if [ ! -s /home/projects/dp_00007/people/hmon/EUostrea/03_datasets/LocalPca/$PREFIX"_"$LG".beagle.gz" ]; then
		echo "Subsetting "$LG
		zcat $BEAGLE | head -n 1 > /home/projects/dp_00007/people/hmon/EUostrea/03_datasets/LocalPca/$PREFIX"_"$LG".beagle"
		zcat $BEAGLE | grep $LG"_" >> /home/projects/dp_00007/people/hmon/EUostrea/03_datasets/LocalPca/$PREFIX"_"$LG".beagle" &
		COUNT=$(( COUNT + 1 ))
		if [ $COUNT == $N_CORE_MAX ]; then
		wait
		COUNT=0
		fi
	else
		echo $LG" was already subsetted"
	fi
done


wait 

COUNT=0
for LG in `cat $LG_LIST`; do
	if [ ! -s /home/projects/dp_00007/people/hmon/EUostrea/03_datasets/LocalPca/$PREFIX"_"$LG".beagle.gz" ]; then
		echo "Gzipping "$LG
		gzip /home/projects/dp_00007/people/hmon/EUostrea/03_datasets/LocalPca/$PREFIX"_"$LG".beagle" &
		COUNT=$(( COUNT + 1 ))
		if [ $COUNT == $N_CORE_MAX ]; then
		wait
		COUNT=0
		fi
	fi
done

🤝
## Split beagle files into smaller windows, each containing a header and the desired number of SNPs
COUNT=0
for LG in `cat $LG_LIST`; do
	echo "Splitting "$LG
	zcat $BEAGLEDIR$PREFIX"_"$LG".beagle.gz" | tail -n +2 | split -d --lines $SNP - --filter='bash -c "{ zcat ${FILE%.*} | head -n1; cat; } > $FILE"' $BEAGLEDIR$PREFIX"_"$LG".beagle.x" &
	COUNT=$(( COUNT + 1 ))
  if [ $COUNT == $N_CORE_MAX ]; then
	  wait
	  COUNT=0
	fi
done

wait

## Gzip these beagle files
COUNT=0
for LG in `cat $LG_LIST`; do
	echo "Zipping "$LG
	gzip $BEAGLEDIR$PREFIX"_"$LG".beagle.x"* &
	COUNT=$(( COUNT + 1 ))
  if [ $COUNT == $N_CORE_MAX ]; then
	  wait
	  COUNT=0
	fi
done

wait


🤝
## Run pcangsd and prepare the local_pca input. 
cd /home/projects/dp_00007/people/hmon/EUostrea
COUNT=0
for LG in `cat $LG_LIST`; do
	if [ -f $BEAGLEDIR"snp_position14feb23_"$SNP"snp_"$LG".tsv" ]; then
		rm $BEAGLEDIR"snp_position14feb23_"$SNP"snp_"$LG".tsv"
	fi
	if [ -f $BEAGLEDIR"pca_summary14feb23_"$SNP"snp_"$LG".tsv" ]; then
		rm $BEAGLEDIR"pca_summary14feb23_"$SNP"snp_"$LG".tsv"
	fi
	bash $LOCAL_PCA_1 $BEAGLEDIR $PREFIX $LG $PC $SNP $PYTHON $PCANGSD $LOCAL_PCA_2 &
	COUNT=$(( COUNT + 1 ))
  if [ $COUNT == $N_CORE_MAX ]; then
	  wait
	  COUNT=0
	fi
done

🤝
```

#on c2 generate dist
```R
#### Run pc_dist and assemble the output #### 
install.packages("data.table")
devtools::install_github("petrelharp/local_pca/lostruct")
library(lostruct)
library(tidyverse)
# Read the input
pca_summary <- read_tsv("~/Desktop/Scripts/Data/LocalPCA_EUostrea/pca_summary13feb23_1000snp_4pc.14feb23.tsv", col_names = F) 
# Run pc_dist with pca_summary
pca_summary <- as.matrix(pca_summary)
attr(pca_summary, 'npc') <- 2
dist <- pc_dist(pca_summary) #Is taking a long ass time.
write_tsv(as.data.frame(dist), '"03_datasets/LocalPca/run_pc_dist_1000snp_2pc_c2.14feb23.tsv', col_names = F)

attr(pca_summary, 'npc') <- 4
dist2 <- pc_dist(pca_summary) #Is taking a long ass time.
write_tsv(as.data.frame(dist2), '"03_datasets/LocalPca/run_pc_dist_1000snp_4pc_c2.14feb23.tsv', col_names = F)


```
