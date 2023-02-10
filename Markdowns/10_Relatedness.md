
NGSrelate within each location
===========

>We first generated a file with allele frequencies (-doMaf 1) and a file with genotype likelihoods (-GL 1: using Samtools) (-doGlf 3) for all samples (n=582) passing sample filtering and sites filtering criteria using ANGSD, using the **LDpruned SNP list**, with additional filtering of -minIND two third (of the individual in the pop) -hwe_pval 1e-6 -minmaf 0.05 -minMapQ 20 -minQ 20 -SNP_pval 1e-6 **Per Population**. The frequency column from the allele frequency file was then extracted and used together with the genotype likelihood file by NgsRelate

- [NGSrelate within each location](#ngsrelate-within-each-location)
  - [Allele frequencies (mafs.gz) and genotype likelihoods (glf.gz) files per pop.](#allele-frequencies-mafsgz-and-genotype-likelihoods-glfgz-files-per-pop)
  - [Formatting mafs files for NGSrelate](#formatting-mafs-files-for-ngsrelate)
  - [run NgsRelate](#run-ngsrelate)


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
OUTPUTFOLDER=/home/projects/dp_00007/people/hmon/EUostrea/03_datasets/Relatedness
BASEDIR=/home/projects/dp_00007/people/hmon/EUostrea
REF=/home/projects/dp_00007/people/hmon/AngsdPopStruct/01_infofiles/fileOegenome10scaffoldC3G.fasta

```

## Allele frequencies (mafs.gz) and genotype likelihoods (glf.gz) files per pop.
```bash
POP=("MOLU" "ZECE" "CRES" "ORIS" "CORS" "PONT"  "RIAE" "MORL" "USAM" "TOLL" "COLN" "BARR" "TRAL" "CLEW" "RYAN" "GREV" "WADD" "NISS" "LOGS" "VENO" "HALS" "THIS" "KALV" "HYPP" "LANG" "BUNN" "DOLV" "HAUG" "HAFR" "INNE" "VAGS" "AGAB" "OSTR")

POP=("HAFR")

for query in ${POP[*]}
do
    N_IND=`cat $BASEDIR/01_infofiles/Jan23--EUostrea_${query}-Fst.list | wc -l`
    angsd -nThreads 40 \
    -bam $BASEDIR/01_infofiles/Jan23--EUostrea_${query}-Fst.list \
    -ref $REF \
    -GL 1 -doGlf 3 -doMaf 1 -minmaf 0.05 -doHWE 1 -doMajorMinor 3 \
    -sites $SNP_LIST -rf $LG_LIST \
    -minMapQ 20 -minQ 20 \
    -minInd $((N_IND*2/3)) -dosnpstat 1 -hwe_pval 1e-6 -SNP_pval 1e-6 \
    -out $OUTPUTFOLDER/3feb23_ngsrelate_minind0.75_LDprunedList_${query}
done
```
🤝
## Formatting mafs files for NGSrelate
>Extract the frequency column from the allele frequency file and remove the header (to make it in the format NgsRelate needs)

```bash
for query in ${POP[*]}
do
  zcat $OUTPUTFOLDER/3feb23_ngsrelate_minind0.75_LDprunedList_${query}.mafs.gz | cut -f6 | sed 1d > $OUTPUTFOLDER/3feb23_ngsrelate_minind0.75_LDprunedList_${query}.freq
done
```
## run NgsRelate

```bash
for query in ${POP[*]}
do
  N_IND=`cat $BASEDIR/01_infofiles/Jan23--EUostrea_${query}-Fst.list | wc -l`
  $NGSRELATE \
  -g $OUTPUTFOLDER/3feb23_ngsrelate_minind0.75_LDprunedList_${query}.glf.gz \
  -f $OUTPUTFOLDER/3feb23_ngsrelate_minind0.75_LDprunedList_${query}.freq \
  -n $((N_IND*1)) -O 3feb23_ngsrelate_minind0.75_LDprunedList_${query}.res \
  -p 10
done
```