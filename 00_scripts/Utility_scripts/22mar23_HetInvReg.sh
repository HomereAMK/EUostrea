





# Runs _ANGSD with depth filters to get the mafs file
N_IND=`cat /home/projects/dp_00007/people/hmon/EUostrea/01_infofiles/bamlist_EUostrea.txt | wc -l`
REF=/home/projects/dp_00007/people/hmon/AngsdPopStruct/01_infofiles/fileOegenome10scaffoldC3G.fasta
BAMLIST=/home/projects/dp_00007/people/hmon/EUostrea/01_infofiles/bamlist_EUostrea.txt

REG=("Reg04" "Reg05" "Reg08")

for query in ${REG[*]}

do
angsd -nThreads 40 -ref $REF -bam $BAMLIST -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
-doCounts 1 -dumpCounts 2 \
-minMapQ 20 -minQ 20 -setMinDepthInd 1 -setMinDepth 600 -setMaxDepth 1200 \
-GL 1 -doMajorMinor 1 -doMaf 1 \
-rf /home/projects/dp_00007/people/hmon/EUostrea/01_infofiles/${query}.txt \
-out /home/projects/dp_00007/people/hmon/EUostrea/03_datasets/Het/Mar23--Inv${query}_NICdepthfilters_het_angsd
done

for query in Reg04 Reg04 Reg04
do
