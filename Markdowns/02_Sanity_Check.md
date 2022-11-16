>SanityCheck with replicates and all the bams, require a new angsd job that will run only on the first chromosome with the replicates
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

    angsd \
     -b $BAMLIST -ref $REF -out $OUTPUTDIR/SC_chr1_Tyler_minMapQ20_minInd0.25_setMinDepthInd1_setMinDepth7_setMaxDepth20000_SNPpval1e-6_minMaf0.05_nov22 \
    -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 \
    -minMapQ 20 -minQ 20 -minInd $((N_IND*1/4)) -setMinDepthInd 1 -setMinDepth 7 -setMaxDepth 20000 \
    -doCounts 1 -dumpCounts 2 \
    -GL 1 -doGlf 2 \
    -doMajorMinor 1 -doMaf 1 -SNP_pval 1e-6 -minMaf 0.05 -rmTriallelic 0.05 -doPost 1 -doGeno 3 -postCutoff 0.95 \
    -doCov 1 -makeMatrix 1 \
    -nThreads 40
                                                                                                                                                                                                                                                                                       