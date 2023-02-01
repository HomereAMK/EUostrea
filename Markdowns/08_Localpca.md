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
#variables
BEAGLE=/home/projects/dp_00007/people/hmon/EUostrea/03_datasets/PopulationStructure/30jan23_prunedLDminweight0.5_PopStruct.beagle.gz
LG_LIST=/home/projects/dp_00007/people/hmon/EUostrea/01_infofiles/List_scaffold_28jan23.txt
N_CORE_MAX=40 # Maximum number of threads to use simulatenously
OUTPUTFOLDER=/home/projects/dp_00007/people/hmon/EUostrea/03_datasets/LocalPca
PREFIX=`echo $BEAGLE | sed 's/\..*//' | sed -e 's#.*/\(\)#\1#'`
BEAGLEDIR=`echo $BEAGLE | sed 's:/[^/]*$::' | awk '$1=$1"/"'`
```

## Subset the beagle file: subset_beagle_by_lg.sh
# This script is used to subset a genome-wide beagle file into smaller files by linkage groups or chromosomes. 
```bash
BEAGLE=/home/projects/dp_00007/people/hmon/EUostrea/03_datasets/PopulationStructure/30jan23_prunedLDminweight0.5_PopStruct.beagle.gz
LG_LIST=/home/projects/dp_00007/people/hmon/EUostrea/01_infofiles/List_scaffold_28jan23.txt
N_CORE_MAX=40 # Maximum number of threads to use simulatenously
OUTPUTFOLDER=/home/projects/dp_00007/people/hmon/EUostrea/03_datasets/LocalPca
PREFIX=`echo $BEAGLE | sed 's/\..*//' | sed -e 's#.*/\(\)#\1#'`
BEAGLEDIR=`echo $BEAGLE | sed 's:/[^/]*$::' | awk '$1=$1"/"'`

COUNT=0
for LG in `cat $LGLIST`; do
	if [ ! -s $OUTPUTFOLDER/$PREFIX"_"$LG".beagle.gz" ]; then
		echo "Subsetting "$LG
		zcat $BEAGLE | head -n 1 > $OUTPUTFOLDER/$PREFIX"_"$LG".beagle"
		zcat $BEAGLE | grep $LG"_" >> $OUTPUTFOLDER/$PREFIX"_"$LG".beagle" &
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
for LG in `cat $LGLIST`; do
	if [ ! -s $OUTPUTFOLDER/$PREFIX"_"$LG".beagle.gz" ]; then
		echo "Gzipping "$LG
		gzip $OUTPUTFOLDER/$PREFIX"_"$LG".beagle" &
		COUNT=$(( COUNT + 1 ))
		if [ $COUNT == $N_CORE_MAX ]; then
		wait
		COUNT=0
		fi
	fi
done
```

## window of 1000snps
```bash
#for 1000snpsWindow
SNP=1000 ## Number of SNPs to include in each window
PC=2 ## Number of PCs to keep for each window
N_CORE_MAX=40 # Maximum number of threads to use simulatenously
COUNT=0
for LG in `cat $LGLIST`; do
	echo "Splitting "$LG
	zcat $PREFIX"_"$LG".beagle.gz" | tail -n +2 | split -d --lines $SNP - --filter='bash -c "{ zcat ${FILE%.*} | head -n1; cat; } > $FILE"' $PREFIX"_"$LG".beagle.x" &
	COUNT=$(( COUNT + 1 ))
if [ $COUNT == $N_CORE_MAX ]; then
	wait
	COUNT=0
	fi
done
````

## Gzip these beagle files
COUNT=0
for LG in `cat $LGLIST`; do
	echo "Zipping "$LG
	gzip $PREFIX"_"$LG".beagle.x"* &
	COUNT=$(( COUNT + 1 ))
  if [ $COUNT == $N_CORE_MAX ]; then
	  wait
	  COUNT=0
	fi
done
