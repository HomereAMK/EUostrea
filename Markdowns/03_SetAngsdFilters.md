Best Practices to Set up the Ansgd filters for Variant calling and more
(This markdown has been written by R. Nicolas Lou with minor modifications from me)
================
- [(This markdown has been written by R. Nicolas Lou with minor modifications from me)](#this-markdown-has-been-written-by-r-nicolas-lou-with-minor-modifications-from-me)
  - [Genome statistics](#genome-statistics)
  - [Establish SNP calling filters](#establish-snp-calling-filters)
      - [Count read depth per position](#count-read-depth-per-position)
  - [SNP variant calling](#snp-variant-calling)
- [Trial dataset](#trial-dataset)
- [Shuffle the IDs in the file](#shuffle-the-ids-in-the-file)
- [Select the first 30 shuffled IDs](#select-the-first-30-shuffled-ids)
- [Run Angsd with the trial](#run-angsd-with-the-trial)


``` r
library(tidyverse)
library(cowplot)
library(RcppCNPy)
library(scales)
library(RColorBrewer)
library(knitr)
library(ade4)
```

## Genome statistics 
Loads modules:
```
module load samtools/1.14
module load gcc/11.1.0
module load intel/perflibs/2020_update4E
module load R/4.2.0
module load ngsadmix/32
```
Counts the number of base pairs in a genome:
```
REF=/home/projects/dp_00007/people/hmon/AngsdPopStruct/01_infofiles/fileOegenome10scaffoldC3G.fasta
samtools faidx $REF
```
Gets a CHRs-Scaffolds list:
```
cat $REF.fai | awk '{sum+=($2)} END {print "Total: " sum}'
grep "scaffold" $REF.fai | awk '{sum+=($2)} END {print "CHRs: " sum}'
grep -v "scaffold" $REF.fai | awk '{sum+=($2)} END {print "Scaffolds: " sum}'
```
| --- | # of Scaffolds | # of bases  | % of genome |
| :---: | :---: | :---: | :---: |
| CHRs | 10 | 858073199 | 100% |

## Establish SNP calling filters

#### Count read depth per position
Here we are interested in obtaining a read depth distribution per position obtained from ANGSD to establish depth filters for SNP calling.
We use `ANGSD -doDepth 1`. More filters can be used with this method, including the minimum base quality score (`-minQ 20`),`-remove_bads 1`, `-only_proper_pairs 1`. Another advantage of this approach is that a list of sites that have \>0 depth across all individuals can be generated (through `-doMAF 1`), and this list can be further filtered through R (based on depth and number of individuals) and can be used for heterozygosity estimation (since non-variable sites are not removed)!

``` bash
#angsd
module load tools computerome_utils/2.0
#module load htslib/1.9
module load bedtools/2.30.0
module load pigz/2.3.4
module load parallel/20210722
#module load angsd/0.929
#variables
REF=/home/projects/dp_00007/people/hmon/AngsdPopStruct/01_infofiles/fileOegenome10scaffoldC3G.fasta
BAMLIST=/home/projects/dp_00007/people/hmon/EUostrea/01_infofiles/bamlist_EUostrea.txt
BASEDIR=/home/projects/dp_00007/people/hmon/EUostrea
OUTPUTFOLDER=/home/projects/dp_00007/people/hmon/EUostrea/03_datasets/SetAngsdFilters

angsd \
-bam $BAMLIST \
-ref $REF \
-out $OUTPUTFOLDER/Nlou_filtered_minq20_minmapq20_angsd0.929_htslib1.9 \
-GL 1 -doMajorMinor 1 -doMaf 1 \
-doDepth 1 -doCounts 1 -maxDepth 6000 -dumpCounts 1 \
-minMapQ 20 -minQ 20 \
-remove_bads 1 -only_proper_pairs 1 \
-nThreads 40
```


## SNP variant calling
```
#angsd
module load tools computerome_utils/2.0
module load htslib/1.16
module load bedtools/2.30.0
module load pigz/2.3.4
module load parallel/20210722
module load angsd/0.940

#variables
REF=/home/projects/dp_00007/people/hmon/AngsdPopStruct/01_infofiles/fileOegenome10scaffoldC3G.fasta
BAMLIST=/home/projects/dp_00007/people/hmon/EUostrea/01_infofiles/bamlist_EUostrea.txt
BASEDIR=/home/projects/dp_00007/people/hmon/EUostrea
OUTPUTFOLDER=/home/projects/dp_00007/people/hmon/EUostrea/03_datasets/VariantCalling
N_IND=`cat /home/projects/dp_00007/people/hmon/EUostrea/01_infofiles/bamlist_EUostrea.txt | wc -l`

```

```
angsd \
-bam $BAMLIST \
-ref $REF \
-out $OUTPUTFOLDER/angsd0.937_htslib1.16_minMapQ20minQ20_minInd145.25_setMinDepthInd5_setMinDepth600setMaxDepth1200 \
-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
-minMapQ 20 -minQ 20 -minInd $((N_IND*1/4)) -setMinDepthInd 5 -setMinDepth 600 -setMaxDepth 1200 \
-doCounts 1 -dumpCounts 2 \
-GL 1 -doGlf 2 \
-doMajorMinor 1 -doMaf 1 -SNP_pval 1e-6 -minMaf 0.05 -rmTriallelic 0.05 -doPost 1 -doGeno 8 \
-doIBS 1 -doCov 1 -makeMatrix 1 \
-nThreads 40
````
        -> Total number of sites analyzed: 840354420
        -> Number of sites retained after filtering: 32477
        [ALL done] cpu-time used =  698796.44 sec
        [ALL done] walltime used =  568548.00 sec

```
```
#Get the annotation file 
cat /home/projects/dp_00007/people/hmon/EUostrea/01_infofiles/bamlist_EUostrea.2.labels | awk '{split($0,a,"_"); print $1"\t"a[1]}' > /home/projects/dp_00007/people/hmon/EUostrea/01_infofiles/bamlist_EUostrea.annot
```

# Trial dataset
# Shuffle the IDs in the file
  sort -R 01_infofiles/bamlist_EUostrea.txt > 01_infofiles/bamlist_EUostreashuffled_ids.txt
# Select the first 30 shuffled IDs
  head -n 30 01_infofiles/bamlist_EUostreashuffled_ids.txt > 01_infofiles/Trialselected_ids_EUostrea.txt
# Run Angsd with the trial
#variables
  REF=/home/projects/dp_00007/people/hmon/AngsdPopStruct/01_infofiles/fileOegenome10scaffoldC3G.fasta
  BAMLIST2=/home/projects/dp_00007/people/hmon/EUostrea/01_infofiles/Trialselected_ids_EUostrea.txt
  OUTPUTFOLDER2=/home/projects/dp_00007/people/hmon/EUostrea/03_datasets/Trial
  N_IND=`cat /home/projects/dp_00007/people/hmon/EUostrea/01_infofiles/bamlist_EUostrea.txt | wc -l`

  angsd \
  -bam $BAMLIST2 \
  -ref $REF \
  -out $OUTPUTFOLDER2/Trial0.940_minMapQ20minQ20_minInd25perc_setMinDepthInd5_setMinDepth600setMaxDepth1200 \
  -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
  -minMapQ 20 -minQ 20 -minInd $((N_IND*1/4)) -setMinDepthInd 5 -setMinDepth 600 -setMaxDepth 1200 \
  -doCounts 1 -dumpCounts 2 \
  -GL 1 -doGlf 2 \
  -doMajorMinor 1 -doMaf 1 -SNP_pval 1e-6 -minMaf 0.05 -rmTriallelic 0.05 -doPost 1 -doGeno 8 \
  -doIBS 1 -doCov 1 -makeMatrix 1 \
  -nThreads 40
