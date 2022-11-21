>SanityCheck with replicates and all the bams, require a new angsd job that will run only on the 1st chromosome with the replicates
    #angsd
    module load tools computerome_utils/2.0
    module load htslib/1.16
    module load bedtools/2.30.0
    module load pigz/2.3.4
    module load parallel/20210722
    module load angsd/0.937

    REF=/home/projects/dp_00007/people/hmon/AngsdPopStruct/01_infofiles/fileOegenome10scaffoldC3G.fasta
    BAMLIST=/home/projects/dp_00007/people/hmon/EUostrea/01_infofiles/bamlist_EUostrea_replicates.txt
    BASEDIR=/home/projects/dp_00007/people/hmon/EUostrea
    OUTPUTFOLDER=/home/projects/dp_00007/people/hmon/EUostrea/03_datasets/SanityCheck
    N_IND=`cat /home/projects/dp_00007/people/hmon/EUostrea/01_infofiles/bamlist_EUostrea_replicates.txt | wc -l`

Parameter | Meaning |
--- | --- |
-setMinDepthInd 1 | Minimum number of reads to consider an individual as having non-missing data is 1 |
-minInd 25% | use only sites where at least 25% individuals have data |
-setMinDepth 7 | minimum total site depth is 7  |
-setMaxDepth 20.000 | maximum total site depth is 20.000|

    angsd \
    -b $BAMLIST -ref $REF -out $OUTPUTFOLDER/SC_chr1_Tyler_minMapQ20_minInd0.25_setMinDepthInd1_setMinDepth7_rmTriallelic0.05minMaf0.05__setMaxDepth20000_SNPpval1e-6_minMaf0.05_nov22 \
    -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 \
    -minMapQ 20 -minQ 20 -minInd $((N_IND*1/4)) -setMinDepthInd 1 -setMinDepth 7 -setMaxDepth 20000 \
    -doCounts 1 -dumpCounts 2 \
    -GL 1 -doGlf 2 \
    -doMajorMinor 1 -doMaf 1 -SNP_pval 1e-6 -minMaf 0.05 -rmTriallelic 0.05 -doPost 1 -doGeno 3 \
    -doIBS 1 -doCov 1 -makeMatrix 1 \
    -nThreads 40 \
    -r scaffold1

SNPs |  Number of sites with at least one individual with missing data | Number of sites with at least 10 individuals with missing data |
--- | --- | --- |
1159625 | 1159625 | 10654 |

> Missing Data on the Variant Calling Sanity Check chr 1

#Get the label list from the bam list
awk '{split($0,a,"/"); print a[9]}' $BAMLIST | awk '{split($0,b,"_"); print b[1]"_"b[2]}' > /home/projects/dp_00007/people/hmon/Flat_oysters/01_infofiles/Bam_list_13dec21.labels

    zcat SC_chr1_Tyler_minMapQ20_minInd0.25_setMinDepthInd1_setMinDepth7_rmTriallelic0.05minMaf0.05__setMaxDepth20000_SNPpval1e-6_minMaf0.05_nov22.counts.gz | tail -n +2 | gawk ' {for (i=1;i<=NF;i++){a[i]+=$i;++count[i]}} END{ for(i=1;i<=NF;i++){print a[i]/count[i]}}' | paste /home/projects/dp_00007/people/hmon/Flat_oysters/01_infofiles/Bam_list_13dec21.labels - > /home/projects/dp_00007/people/hmon/Flat_oysters/02_ngsLDOutput/Dataset_I/Leona20dec21_SNPs_11jan22.GL-RealCoverage.txt
