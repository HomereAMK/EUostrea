#!/bin/bash
## This script is the first dependency of /workdir/genomic-data-analysis/scripts 
## It will loop through all windowed beagle files in each LG
## For each beagle file, it runs pcangsd first, and then runs an R script (/workdir/genomic-data-analysis/scripts/local_pca_2.R) to process the covariance matrix.  

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

BEAGLEDIR=$1
PREFIX=$2
LG=$3
PC=$4
SNP=$5
PYTHON=$6
PCANGSD=$7
LOCAL_PCA_2=$8

## Set maximum number of threads to 1
export OMP_NUM_THREADS=1

## Loop through each windowed beagle file in the same linkage group (or chromosome)
for INPUT in `ls $BEAGLEDIR$PREFIX"_"$LG".beagle.x"*".gz"`; do
	## Run pcangsd
	pcangsd -beagle $INPUT -o $INPUT --threads 5
	## Process pcangsd output
	Rscript --vanilla $LOCAL_PCA_2 $INPUT".cov" $PC $SNP $INPUT $LG $BEAGLEDIR
done
