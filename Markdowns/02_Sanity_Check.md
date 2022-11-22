## SanityCheck with replicates and all the bams, require a new angsd job that will run only on the 1st chromosome with the replicates 


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


## Missing Data on the Variant Calling Sanity Check chr 1
    #Get the annotation file 
    cat $BASEDIR/01_infofiles/bamlist_EUostrea_replicates.labels  | awk '{split($0,a,"_"); print $1"\t"a[1]}' > $BASEDIR/01_infofiles/bamlist_EUostrea_replicates.annot

    #Gets Real Coverage (Genotype Likelihoods):
    zcat $OUTPUTFOLDER/SC_chr1_Tyler_minMapQ20_minInd0.25_setMinDepthInd1_setMinDepth7_rmTriallelic0.05minMaf0.05__setMaxDepth20000_SNPpval1e-6_minMaf0.05_nov22.counts.gz | tail -n +2 | gawk ' {for (i=1;i<=NF;i++){a[i]+=$i;++count[i]}} END{ for(i=1;i<=NF;i++){print a[i]/count[i]}}' | paste $BASEDIR/01_infofiles/bamlist_EUostrea_replicates.labels - > $OUTPUTFOLDER/SC_chr1_Tyler_minMapQ20_minInd0.25_setMinDepthInd1_setMinDepth7_rmTriallelic0.05minMaf0.05__setMaxDepth20000_SNPpval1e-6_minMaf0.05_nov22.GL-RealCoverage.txt

    #Gets Missing Data (Genotype Likelihoods):
    N_SITES=`zcat /home/projects/dp_00007/people/hmon/EUostrea/03_datasets/SanityCheck/SC_chr1_Tyler_minMapQ20_minInd0.25_setMinDepthInd1_setMinDepth7_rmTriallelic0.05minMaf0.05__setMaxDepth20000_SNPpval1e-6_minMaf0.05_nov22.beagle.gz | tail -n +2 | wc -l`

    zcat $OUTPUTFOLDER/SC_chr1_Tyler_minMapQ20_minInd0.25_setMinDepthInd1_setMinDepth7_rmTriallelic0.05minMaf0.05__setMaxDepth20000_SNPpval1e-6_minMaf0.05_nov22.beagle.gz | tail -n +2 | perl /home/projects/dp_00007/apps/Scripts/call_geno.pl --skip 3 | cut -f 4- | awk '{ for(i=1;i<=NF; i++){ if($i==-1)x[i]++} } END{ for(i=1;i<=NF; i++) print i"\t"x[i] }' | paste $BASEDIR/01_infofiles/bamlist_EUostrea_replicates.labels - | awk -v N_SITESawk="$N_SITES" '{print $1"\t"$3"\t"$3*100/N_SITESawk}' > /$OUTPUTFOLDER/SC_chr1_Tyler_minMapQ20_minInd0.25_setMinDepthInd1_setMinDepth7_rmTriallelic0.05minMaf0.05__setMaxDepth20000_SNPpval1e-6_minMaf0.05_nov22.GL-MissingData.txt

![](https://github.com/HomereAMK/EUostrea/blob/main/Figures/SanityCheck/DendroIbsMat_SC_chr1_Tyler_minMapQ20_minInd0.25_setMinDepthInd1_setMinDepth7_rmTriallelic0.05minMaf0.05__setMaxDepth20000_SNPpval1e-6_minMaf0.05_nov22.pdf)<!-- -->

