
NGSrelate within each location
===========

>We first generated a file with allele frequencies (-doMaf 1) and a file with genotype likelihoods (-doGlf 3) for all samples (n=582) passing sample filtering and sites filtering criteria using ANGSD, using the **LDpruned SNP list**, with additional filtering of -hwe_pval 1e-6 -minmaf 0.05 -minMapQ 20 -minQ 20 -SNP_pval 1e-6 **Per Population**. The frequency column from the allele frequency file was then extracted and used together with the genotype likelihood file by NgsRelate

- [NGSrelate within each location](#ngsrelate-within-each-location)
  - [Allele frequencies (mafs.gz) and genotype likelihoods (glf.gz) files per pop.](#allele-frequencies-mafsgz-and-genotype-likelihoods-glfgz-files-per-pop)


```bash
#angsd
module load tools computerome_utils/2.0
module load htslib/1.16
module load bedtools/2.30.0
module load pigz/2.3.4
module load parallel/20210722
module load angsd/0.940

#variables
NGSRELATE=/home/projects/dp_00007/people/hmon/Software/ngsRelate/ngsRelate
SNP_LIST=/home/projects/dp_00007/people/hmon/EUostrea/03_datasets/LDpruning/LDprunedlist_rightmafs_AllCHRs.min_weight0.5_23jan23
LG_LIST=/home/projects/dp_00007/people/hmon/EUostrea/01_infofiles/List_scaffold_28jan23.txt
```

## Allele frequencies (mafs.gz) and genotype likelihoods (glf.gz) files per pop.
