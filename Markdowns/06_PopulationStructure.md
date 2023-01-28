Population Structure Analysis
================

  - [Modules](#Modules)
  - [Re-run angsd LD Pruned SNPs list minweight0.5](#re-run-angsd-LD-Pruned-snps-list-minweight05)


## Modules
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

## Re-run angsd LD Pruned SNPs list minweight0.5
REF=/home/projects/dp_00007/people/hmon/AngsdPopStruct/01_infofiles/fileOegenome10scaffoldC3G.fasta
BAMLIST=/home/projects/dp_00007/people/hmon/EUostrea/01_infofiles/bamlist_EUostrea.txt
BASEDIR=/home/projects/dp_00007/people/hmon/EUostrea
N_IND=`cat /home/projects/dp_00007/people/hmon/EUostrea/01_infofiles/bamlist_EUostrea.txt | wc -l`
OUTPUTFOLDER=/home/projects/dp_00007/people/hmon/EUostrea/03_datasets/PopulationStructure

  angsd \
  -bam $BAMLIST \
  -ref $REF \
  -sites \
  -out $OUTPUTFOLDER \
  -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
  -minMapQ 20 -minQ 20 -setMinDepthInd 1 -setMinDepth 600 -setMaxDepth 1200 \
  -doCounts 1 -dumpCounts 2 \
  -GL 1 -doGlf 2 \
  -doMajorMinor 1 -doMaf 1 -SNP_pval 1e-6 -minMaf 0.05 -rmTriallelic 0.05 -doPost 1 -doGeno 8 \
  -doIBS 1 -doCov 1 -makeMatrix 1 \
  -nThreads 40 \
  -P 40
