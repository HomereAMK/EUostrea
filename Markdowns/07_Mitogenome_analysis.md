# MitOyster

Recent evolutionary history of Scandinavian flat oysters is unclear. Both phylogenetic relationship and demographic history remain fuzzy given the extensive historical and present translocations this commercially valuable species is the object of.
In a recent publications, Hayer and colleagues provided historical haplotypes coming from 19th century _Ostrea edulis __ ____________
mitochondrial genomes as well as a new mitondrial genome refence coming from a present (2018) Limfjorden (Denmark) flat oyster.

mtDNA has a much higher mutation rate and shorter coalescence times than most nuclear regions, valuable for the study of more recent demographic processes at an intraspecific level (Brown et al. 1979; Avise et al. 1987). mitochondrial genome sequences extracted from shotgun HTS can provide a valuable perspective into the phylogeography of a species that is complementary to the nuclear signals, and can form a link between the past studies based on short mitochondrial sequences and future studies using whole nuclear genomes.

MtGenome reference: Hayer et al., 2021 MT663266 (Limfjorden) 16,356bp
-> no need to go prepare a new assembly 

I. Processing reads

## I. Mapping with BWA

[II Analyses](II-Analyses)

## II Analyses
Regional structure and phylogenetics
>Angsd to get a beagle file 
```
# Load module angsd
# Load module phylogeny
# Load module perl
module load tools
module load ngs
module load htslib/1.16
module load angsd/0.940
module load gcc/8.2.0
module load openmpi/gcc/64
module load raxml-ng/1.1.0
module load ngsdist/20191118
module load fastme/2.1.6.2
module load raxml/8.2.11
module load perl/5.30.2
module load samtools/1.16
```


```
#Consensus Fasta seq for a all pop or 
POP_BAMLIST=/home/projects/dp_00007/people/hmon/MitOyster/01_infofiles/List_phylogenyMT_7jun22.txt
N_IND=`cat /home/projects/dp_00007/people/hmon/MitOyster/01_infofiles/List_phylogenyMT_7jun22.txt | wc -l`
OUTPUTFOLDER=/home/projects/dp_00007/people/hmon/MitOyster/03_results/scandinavia
OUTNAME=30jan23_UltraMT
GENOME="01_infofiles/MT663266.fasta"

angsd -bam $POP_BAMLIST -ref $GENOME -out $OUTPUTFOLDER/$OUTNAME \
-minInd $((N_IND*2/3)) \
-remove_bads 1 -uniqueOnly 1 -baq 1 -C 50 -minMapQ 30 -minQ 30 \
-doCounts 1 -GL 1 -doGlf 2 -doMajorMinor 1 -doMaf 1 -doPost 2 -doGeno 3 -doPlink 2 -geno_minDepth 3 -setMaxDepth $((323*323)) -dumpCounts 2 -postCutoff 0.95 -doHaploCall 1 -doBcf 1 
```

>phylogeny
```bash
awk '{split($0,a,"/"); print a[10]}' $POP_BAMLIST | awk '{split($0,b,"_"); print b[1]"_"b[2]}' > /home/projects/dp_00007/people/hmon/MitOyster/01_infofiles/List_phylogenyMT_7jun22.labels
```

>31jan23 quick phylogeny 
```bash
zcat $OUTPUTFOLDER/$OUTNAME.haplo.gz | cut -f 4- | tail -n +2 | perl /home/projects/dp_00007/people/hmon/MitOyster/00_scripts/tsv_merge.pl --transp --ofs '' - | awk 'NR==FNR{id=$1; sub(".*\\/","",id); sub("\\..*","",id); x[FNR]=id} NR!=FNR{ print ">"x[FNR]"\n"$1}' /home/projects/dp_00007/people/hmon/MitOyster/01_infofiles/List_phylogenyMT_7jun22.labels - > $OUTPUTFOLDER/$OUTNAME.fasta
```
>Raxml job
```bash
raxml-ng --threads 10 --search --model GTR+F --site-repeats on --msa $OUTPUTFOLDER/$OUTNAME.fasta --prefix $OUTPUTFOLDER/$OUTNAME.Possible.raxmlng_31jan23
```

First Tree done on https://itol.embl.de/tree/192389368159681675247682 