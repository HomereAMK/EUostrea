Best Practices to Set up the Ansgd filters for Variant calling and more
(This markdown has been written by Nicolas Lou with minor modifications from me)
================

-   [Establish SNP calling filters](#establish-snp-calling-filters)
-   [Reference bias filtering](#reference-bias-filtering)
-   [SNP calling](#snp-calling)
-   [Generate a bamlist for each population](#generate-a-bamlist-for-each-population)
-   [Get minor allele frequencies per population](#get-minor-allele-frequencies-per-population)
-   [Estimate Fst in each pair of populations](#estimate-fst-in-each-pair-of-populations)
-   [Get “effective sample size” per sample per site](#get-effective-sample-size-per-sample-per-site)


``` r
library(tidyverse)
library(cowplot)
library(RcppCNPy)
library(scales)
library(RColorBrewer)
library(knitr)
library(ade4)
```

## Establish SNP calling filters

#### Count read depth per position
Here we are interested in obtaining a read depth distribution per position obtained from ANGSD to establish depth filters for SNP calling.
We use `ANGSD -doDepth 1`. More filters can be used with this method, including the minimum base quality score (`-minQ 20`),`-remove_bads 1`, `-only_proper_pairs 1`. Another advantage of this approach is that a list of sites that have \>0 depth across all individuals can be generated (through `-doMAF 1`), and this list can be further filtered through R (based on depth and number of individuals) and can be used for heterozygosity estimation (since non-variable sites are not removed)!

``` bash
#angsd
module load tools computerome_utils/2.0
module load htslib/1.16
module load bedtools/2.30.0
module load pigz/2.3.4
module load parallel/20210722
module load angsd/0.937
#variables
REF=/home/projects/dp_00007/people/hmon/AngsdPopStruct/01_infofiles/fileOegenome10scaffoldC3G.fasta
BAMLIST=/home/projects/dp_00007/people/hmon/EUostrea/01_infofiles/bamlist_EUostrea.txt
BASEDIR=/home/projects/dp_00007/people/hmon/EUostrea
OUTPUTFOLDER=/home/projects/dp_00007/people/hmon/EUostrea/03_datasets/SetAngsdFilters

angsd \
-bam $BAMLIST \
-ref $REF \
-out $OUTPUTFOLDER/Nlou_filtered_minq20_minmapq20 \
-GL 1 -doMajorMinor 1 -doMaf 1 \
-doDepth 1 -doCounts 1 -maxDepth 2000 -dumpCounts 1 \
-minMapQ 20 -minQ 20 \
-remove_bads 1 -only_proper_pairs 1 \
-nThreads 40
```


``` r
depth_hist_angsd <- read_lines("../../Nlou_filtered_minq20_minmapq20.depthGlobal") %>%
  str_split(pattern = "\t") %>%
  .[[1]] %>%
  .[1:2001] %>%
  as.integer() %>%
  tibble(depth=0:2000, count=., method="angsd")

histo_distrib_<- depth_hist_angsd %>%
  filter(depth>0, depth<1000) %>%
  ggplot(aes(x=depth, y=count, color=method)) +
  geom_freqpoly(stat = "identity") +
  theme_cowplot()
```
