

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
#-finished

# Generates a `.bed` file based on the `.mafs` file:
for query in ${REG[*]}
do
zcat  /home/projects/dp_00007/people/hmon/EUostrea/03_datasets/Het/Mar23--Inv${query}_NICdepthfilters_het_angsd.mafs.gz | cut -f1,2 | tail -n +2 | awk '{print $1"\t"$2-1"\t"$2}' | bedtools merge -i - >  /home/projects/dp_00007/people/hmon/EUostrea/03_datasets/Het/Mar23--Inv${query}_NICdepthfilters_het_angsd.bed
done

#-finished

# Creates a position file based on this new `.bed`:
for query in ${REG[*]}
do
awk '{print $1"\t"($2+1)"\t"$3}'  /home/projects/dp_00007/people/hmon/EUostrea/03_datasets/Het/Mar23--Inv${query}_NICdepthfilters_het_angsd.bed > /home/projects/dp_00007/people/hmon/EUostrea/03_datasets/Het/Mar23--Inv${query}_NICdepthfilters_het_angsd.bed.pos
done

#-finished

# Indexs the `.pos` file created above:
for query in ${REG[*]}
do
angsd sites index /home/projects/dp_00007/people/hmon/EUostrea/03_datasets/Het/Mar23--Inv${query}_NICdepthfilters_het_angsd.bed.pos
done

#-finished

# Runs _ANGSD_ under -doSaf::
for query in ${REG[*]}
do
parallel --plus angsd -i {} -ref $REF -anc $REF -sites /home/projects/dp_00007/people/hmon/EUostrea/03_datasets/Het/Mar23--Inv${query}_NICdepthfilters_het_angsd.bed.pos -rf /home/projects/dp_00007/people/hmon/EUostrea/01_infofiles/${query}.txt -GL 1 -doSaf 1 -remove_bads 1 -uniqueOnly 1 -baq 1 -C 50 -minMapQ 20 -minQ 20 -out /home/projects/dp_00007/data/hmon/angsd_Het/${query}/${query}{/...} :::: /home/projects/dp_00007/people/hmon/EUostrea/01_infofiles/bamlist_EUostrea.txt
done

#-finished


# Calculates fractions:
for query in ${REG[*]}
do
parallel --tmpdir /home/projects/dp_00007/data/hmon/angsd_Het/${query} -j 1 --plus "realSFS -fold 1 -P 10 {} > /home/projects/dp_00007/data/hmon/angsd_Het/${query}/{/..}.het" ::: /home/projects/dp_00007/data/hmon/angsd_Het/${query}/*.saf.idx
done

#-runnning

# Run 40 instances of realSFS in parallel using specified tmpdir and plus options 
for query in ${REG[*]}
do
parallel --tmpdir /home/projects/dp_00007/data/hmon/angsd_Het/${query} --plus "realSFS -fold 1 -P 10 {} > /home/projects/dp_00007/data/hmon/angsd_Het/${query}/{/..}.het" -j 1 --memfree 160000 ::: /home/projects/dp_00007/data/hmon/angsd_Het/${query}/*.saf.idx
done

```
🤝
# Calculates heterozygosity:
```bash
for query in ${REG[*]}
do
cd /home/projects/dp_00007/data/hmon/angsd_Het/${query}
fgrep '.' ${query}*.het | tr ":" " " | awk '{print $1"\t"$3/($2+$3)*100}' | sed -r 's/.het//g' | awk '{split($0,a,"_"); print $1"\t"a[1]"\t"$2"\t"$3'} | awk '{split($2,a,"-"); print $1"\t"a[2]"\t"$3'} > /home/projects/dp_00007/people/hmon/EUostrea/03_datasets/Het/${query}_24mar23--AllSamples_setMinDepth600_setMaxDepth1200--HET.txt
done
```


