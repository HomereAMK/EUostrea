DATASET PROCESSING
================

- [Downsize MORL and USAM fastq files.](#downsize-morl-and-usam-fastq-files)
    - [With BBMAP: samplerate=0.2 Randomly output only this fraction of reads; 1 means sampling is disabled.](#with-bbmap-samplerate02-randomly-output-only-this-fraction-of-reads-1-means-sampling-is-disabled)
- [Fastq + Bam stats](#fastq--bam-stats)


## Downsize MORL and USAM fastq files.
    module load tools ngs  
    module load jdk/19 openjdk/19 java/1.8.0-openjdk jre/1.8.0-openjdk 
    module load bbmap/38.90

### With BBMAP: samplerate=0.2 Randomly output only this fraction of reads; 1 means sampling is disabled.

<details>
<summary> FOR USAM population </summary>
```bash

DIRFQ=/home/projects/dp_00007/people/hmon/Novaseq_MLX_USA

for POP in USA
    do
        for IND in `echo -n 1 2 3 4 5 6 7 8 9` 
        do
            for NUM in `echo -n  103 104 105 106 117 125 137 138 139 141 150 151 152 164 173 177 185 186 187`
            do
            reformat.sh \
            in1=$DIRFQ/NS.1445.002.IDT_i7_${NUM}---IDT_i5_${NUM}.${POP}-${IND}_R1.fastq.gz \
            in2=$DIRFQ/NS.1445.002.IDT_i7_${NUM}---IDT_i5_${NUM}.${POP}-${IND}_R2.fastq.gz \
            out1=$DIRFQ/${POP}M_0${IND}_NS.1445.002.IDT_i7_${NUM}---IDT_i5_${NUM}_DS0.2_1.fq.gz \
            out2=$DIRFQ/${POP}M_0${IND}_NS.1445.002.IDT_i7_${NUM}---IDT_i5_${NUM}_DS0.2_2.fq.gz \
            samplerate=0.2
            done
        done
    done

for POP in USA
    do
        for IND in `echo -n 10 11 12 13 14 15 16 17 18 19` 
        do
            for NUM in `echo -n  103 104 105 106 117 125 137 138 139 141 150 151 152 164 173 177 185 186 187`
            do
            reformat.sh \
            in1=$DIRFQ/NS.1445.002.IDT_i7_${NUM}---IDT_i5_${NUM}.${POP}-${IND}_R1.fastq.gz \
            in2=$DIRFQ/NS.1445.002.IDT_i7_${NUM}---IDT_i5_${NUM}.${POP}-${IND}_R2.fastq.gz \
            out1=$DIRFQ/${POP}M_${IND}_NS.1445.002.IDT_i7_${NUM}---IDT_i5_${NUM}_DS0.2_1.fq.gz \
            out2=$DIRFQ/${POP}M_${IND}_NS.1445.002.IDT_i7_${NUM}---IDT_i5_${NUM}_DS0.2_2.fq.gz \
            samplerate=0.2
            done
        done
    done```

</details>

<details>
<summary> FOR MORL population </summary>   
```bash

for POP in MLX
    do
        for IND in `echo -n 1 2 4 6 7 8 9`
        do
            for NUM in `echo -n 101 113 126 127 128 129 140 149 161 162 163 165 174 175 176 188 189`
            do
            reformat.sh \
            in1=$DIRFQ/NS.1445.002.IDT_i7_${NUM}---IDT_i5_${NUM}.${POP}${IND}_R1.fastq.gz \
            in2=$DIRFQ/NS.1445.002.IDT_i7_${NUM}---IDT_i5_${NUM}.${POP}${IND}_R2.fastq.gz \
            out1=$DIRFQ/MORL_0${IND}_NS.1445.002.IDT_i7_${NUM}---IDT_i5_${NUM}_DS0.2_1.fq.gz \
            out2=$DIRFQ/MORL_0${IND}_NS.1445.002.IDT_i7_${NUM}---IDT_i5_${NUM}_DS0.2_2.fq.gz \
            samplerate=0.2
            done
        done
    done

for POP in MLX
    do
        for IND in `echo -n 10 11 12 13 14 15 16 17 18 19` 
        do
            for NUM in `echo -n 101 113 126 127 128 129 140 149 161 162 163 165 174 175 176 188 189`
            do
            reformat.sh \
            in1=$DIRFQ/NS.1445.002.IDT_i7_${NUM}---IDT_i5_${NUM}.${POP}${IND}_R1.fastq.gz \
            in2=$DIRFQ/NS.1445.002.IDT_i7_${NUM}---IDT_i5_${NUM}.${POP}${IND}_R2.fastq.gz \
            out1=$DIRFQ/MORL_${IND}_NS.1445.002.IDT_i7_${NUM}---IDT_i5_${NUM}_DS0.2_1.fq.gz \
            out2=$DIRFQ/MORL_${IND}_NS.1445.002.IDT_i7_${NUM}---IDT_i5_${NUM}_DS0.2_2.fq.gz \
            samplerate=0.2
            done
        done
    done```
</details>


## Processing USAM MORL fastq
<details>
<summary> click for preprocess code </summary>

```bash 
!/bin/bash
#PBS -d /home/projects/dp_00007/people/hmon/Novaseq_MLX_USA
#PBS -W group_list=dp_00007
#PBS -A dp_00007
#PBS -N trimmomatic__BASE__
#PBS -o __BASE__trim.out
#PBS -e __BASE__trim.err
#PBS -l walltime=00:04:00:00
#PBS -l mem=70g
#PBS -l ncpus=5
#PBS -r n


#Go to the directory from where the job was submitted (initial directory is $HOME)
echo Working directory is $PBS_O_WORKDIR
cd $PBS_O_WORKDIR

#Define number of processors
NPROCS=`wc -l < $PBS_NODEFILE`
echo This job has allocated $NPROCS nodes

#Load all required modules for the job
module load tools
module load ngs
module load jre/1.8.0-openjdk
module load trimmomatic/0.38

#Global variables
base=__BASE__
ADAPTERS=home/projects/dp_00007/people/hmon/Shucking/01_infofiles/NexteraPE_NT.fa

java -jar /services/tools/trimmomatic/0.38/trimmomatic-0.38.jar PE \
        -threads 8 \
        -phred33 \
	"$base"_1.fq.gz \
	"$base"_2.fq.gz \
	"$base"_1.paired.fq.gz \
	"$base"_1.single.fq.gz \
	"$base"_2.paired.fq.gz \
	"$base"_2.single.fq.gz \
	ILLUMINACLIP:$ADAPTERS:2:30:10:1:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:40


#Load all required modules for the job
module load gcc/8.2.0
module load tools
module load ngs
module load bwa/0.7.17
module load samtools/1.12

#Global variables
GENOME="/home/projects/dp_00007/people/hmon/Shucking/01_infofiles/fileOegenome10scaffoldC3G.fasta"
NCPU=8
base=__BASE__

#Align reads
    echo "Aligning $base"
    ID=$(echo "@RG\tID:$base\tSM:$base\tPL:Illumina")

  #Align reads 1 step
    bwa mem -t "$NCPU" \
        -R "$ID" \
        "$GENOME" \
        "$base"_1.paired.fq.gz "$base"_2.paired.fq.gz >"$base".sam

        # Create bam file
    echo "Creating bam for $base"

    samtools view -bS -h -q 20 -F 4 \
    "$base".sam >"$base".bam


     echo "Creating sorted bam for $base"
        samtools sort "$base".bam -o "$base".sort.minq20.bam
        samtools index "$base".sort.minq20.bam

   #Clean up
    echo "Removing "$base".sam"
    echo "Removing "$base".bam"

        rm "$base".sam
        rm "$base".bam


#loading modules
module load tools
module load ngs
module load jre/1.8.0
module load picard-tools/2.25.2
module load parallel/20160822
module load java/1.8.0
module load bamutil/1.0.14

#tryout with NO CIGAR on MarkDuplicates
java -jar /services/tools/picard-tools/2.25.2/picard.jar MarkDuplicates \
I="$base".sort.minq20.bam \
O="$base".nocig.dedup.minq20.bam \
M="$base".duprmmetrics.txt \
REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=SILENT

#scripts ClipOverlap with NO CIGAR on MarkDuplicates
/services/tools/bamutil/1.0.14/bam clipOverlap \
--in "$base".nocig.dedup.minq20.bam \
--out "$base".nocig.dedup_clipoverlap.minq20.bam \
--stats


#ressources
module load tools
module load ngs
module load samtools/1.12
module load parallel/20160822
module load java/1.8.0
module load bamutil/1.0.14
module load gatk/3.8-0
module load jre/1.8.0-openjdk
module load picard-tools/2.25.2

#Index bam files
samtools index "$base".nocig.dedup_clipoverlap.minq20.bam

#Create list of potential in-dels nocig
java -jar /services/tools/gatk/3.8-0/GenomeAnalysisTK.jar \
-T RealignerTargetCreator \
-R $GENOME \
-I "$base".nocig.dedup_clipoverlap.minq20.bam  \
-o "$base".all_samples_for_indel_realigner.nocig.minq20.intervals

#Run the indel realigner tool nocig
java -jar /services/tools/gatk/3.8-0/GenomeAnalysisTK.jar \
-T IndelRealigner \
-R $GENOME \
-I "$base".nocig.dedup_clipoverlap.minq20.bam \
-targetIntervals "$base".all_samples_for_indel_realigner.nocig.minq20.intervals \
--consensusDeterminationModel USE_READS  --nWayOut _minq20.nocig.realigned.bam

```

</details>

## Fastq + Bam stats
<details>
<summary> Counts for USAM and MORL downsized </summary>
```bash 
#Module 
module load tools
module load ngs
module load samtools/1.14

#Global variables
base=__BASE__
DIR=/home/projects/dp_00007/people/hmon/Novaseq_MLX_USA
#raw reads
a=`zcat $DIR/"$base"_1.fq.gz  | wc -l | awk '{print $1/4}'` #raw read forward
b=`zcat $DIR/"$base"_2.fq.gz | wc -l | awk '{print $1/4}'` #raw read reverse
echo $(( $a + $b )) > downS_depth/"$base".count_fastq_1.tmp
#raw bases
c=`zcat $DIR/"$base"_1.fq.gz | awk 'NR%4==2' | tr -d "\n" | wc -m` 
d=`zcat $DIR/"$base"_2.fq.gz | awk 'NR%4==2' | tr -d "\n" | wc -m`
echo $(( $c + $d )) > $DIR/downS_depth/"$base".count_fastq_2.tmp

trim bases
e=`zcat $DIR/"$base"_1.paired.fq.gz | awk 'NR%4==2' | tr -d "\n" | wc -m` 
f=`zcat $DIR/"$base"_2.paired.fq.gz | awk 'NR%4==2' | tr -d "\n" | wc -m` 
echo $(( $e + $f )) > $DIR/downS_depth/"$base".count_fastq_3.tmp 

#mapped bases
samtools stats $DIR/"$base".sort.minq20.bam -@ 12 | grep ^SN | cut -f 2- | grep "^bases mapped (cigar)" | cut -f 2 > $DIR/downS_depth/"$base".count_bam_1.tmp
    
#deduplicate mapped bases
samtools stats $DIR/"$base".nocig.dedup_clipoverlap.minq20.bam -@ 12 | grep ^SN | cut -f 2- | grep "^bases mapped (cigar)" | cut -f 2  > $DIR/downS_depth/"$base".count_bam_2.tmp

#realigned around indels mapped bases
samtools stats $DIR/"$base".nocig.dedup_clipoverlap.minq20_minq20.nocig.realigned.bam -@ 12 | grep ^SN | cut -f 2- | grep "^bases mapped (cigar)" | cut -f 2  > $DIR/downS_depth/"$base".count_bam_3.tmp
#population tag    
echo Novaseq_MLX_USA/"$base"_1.fq.gz |awk '{split($0,a,"_"); print a[2]}' | awk '{split($0,a,"/"); print a[2]}' > $DIR/downS_depth/"$base".count_pop_1.tmp

RAWREADS=`cat $DIR/downS_depth/"$base".count_fastq_1.tmp`
RAWBASES=`cat $DIR/downS_depth/"$base".count_fastq_2.tmp`
ADPTERCLIPBASES=`cat $DIR/downS_depth/"$base".count_fastq_3.tmp`
MAPPEDBASES=`cat $DIR/downS_depth/"$base".count_bam_1.tmp`
DEDUPMAPPEDBASES=`cat $DIR/downS_depth/"$base".count_bam_2.tmp`
REALIGNEDMAPPEDBASES=`cat $DIR/downS_depth/"$base".count_bam_3.tmp`
POP=`cat $DIR/downS_depth/"$base".count_pop_1.tmp`

printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" $base $POP $RAWREADS $RAWBASES $ADPTERCLIPBASES $MAPPEDBASES $DEDUPMAPPEDBASES $REALIGNEDMAPPEDBASES >> $DIR/downS_depth/Summary_DS_USAMMORL_lcWGS_14nov22.txt``` 

</details>

<details>
<summary> Depth stats</summary> 
```bash 
WORKDIR=/home/projects/dp_00007/people/hmon/EUostrea
BAMDIR=/home/projects/dp_00007/people/hmon/Bamfile_EUostrea

#Clean session
cd $WORKDIR
rm 00_scripts/Utility_scripts/DEPTH*sh
#launch scripts for c2 screen
cd $BAMDIR
for file in $(ls *.bam |sed -e 's/.nocig.dedup_clipoverlap.minq20_minq20.nocig.realigned.bam//g'|sort -u)  #only the nocig retry
do
cd $WORKDIR
base=$(basename "$file")
toEval="cat 00_scripts/Utility_scripts/Samtools_depth.sh | sed 's/__BASE__/$base/g'"; eval $toEval > $WORKDIR/00_scripts/Utility_scripts/DEPTH_$base.sh
done

#launch scripts for c2 screen USAM and MORL 
cd /home/projects/dp_00007/people/hmon/Novaseq_MLX_USA/
rm DEPTH*sh
for file in $(ls *nocig.dedup_clipoverlap.minq20_minq20.nocig.realigned.bam |sed -e 's/.nocig.dedup_clipoverlap.minq20_minq20.nocig.realigned.bam//g'|sort -u)  #only the nocig retry
do
    base=$(basename "$file")
    toEval="cat Samtools_depth.sh | sed 's/__BASE__/$base/g'"; eval $toEval > DEPTH_$base.sh
done

#Submit jobs
for i in $(ls DEPTH*sh); do qsub $i; done```

</details>

## Depth plots
<details>
<summary> Make list of depth files to run the Rscript by pop</summary>
```bash
WORKDIR=/home/projects/dp_00007/people/hmon/EUostrea
     
for POP in AGAB BARR BUNN CLEW COLN CORS CRES DOLV GREV HAFR HALS HAUG HYPP INNE KALV LANG LOGS MOLU MORL NISS ORIS OSTR PONT RIAE RYAN THIS TOLL TRAL USAM VAGS VENO WADD ZECE
    do
        BAMSDEPTH="$WORKDIR"/02_data/Depth/${POP}*_depth.gz
        ls $BAMSDEPTH > "$WORKDIR"/01_infofiles/list.${POP}.depth
done```

</details>

<details>
<summary> Calculate % genome covered + Depth per pop </summary>   
    # Load module 
    module load tools
    module load ngs

    ## Load modules FOR R 
    module load gsl/2.6
    module load perl/5.20.1
    module load samtools/1.11
    module load imagemagick/7.0.10-13
    module load gdal/2.2.3
    module load geos/3.8.0
    module load jags/4.2.0
    module load hdf5
    module load netcdf
    module load boost/1.74.0
    module load openssl/1.0.0
    module load lapack
    module load udunits/2.2.26
    module load proj/7.0.0
    module load gcc/10.2.0
    module load intel/perflibs/64/2020_update2
    module load R/4.0.0
    
``` r
    #Clean space
    rm(list=ls())
    #
    library(Rserve)
    library(tidyverse)

    #var
    #for AGAB
    #Loop "BARR", "BUNN", "CLEW", "COLN", "CORS", "CRES", "DOLV", "GREV", "HAFR", "HALS", "HAUG", "HYPP", "INNE", "KALV", "LANG", "LOGS", 
    "MOLU", "MORL", "NISS", "ORIS", "OSTR", "PONT", "RIAE", "RYAN", "THIS", "TOLL", "TRAL", "USAM", "VAGS", "VENO", "WADD", "ZECE"))


    basedir <- "/home/projects/dp_00007/people/hmon/EUostrea" # Make sure to edit this to match your $BASEDIR
    bam_list <- read_lines(paste0(basedir, "/01_infofiles/list.WADD.depth"))

        for (i in 1:length(bam_list)){

        bamfile = bam_list[i]
        #Compute depth stats
        depth <- read_tsv(paste0(bamfile), col_names = F)$X1
        mean_depth <- mean(depth)
        sd_depth <- sd(depth)
        mean_depth_nonzero <- mean(depth[depth > 0])
        mean_depth_within2sd <- mean(depth[depth < mean_depth + 2 * sd_depth])
        median <- median(depth)
        presence <- as.logical(depth)
        proportion_of_reference_covered <- mean(presence)

        #Bind stats into dataframe and store sample-specific per base depth and presence data
        if (i==1){
            output <- data.frame(bamfile, mean_depth, sd_depth, mean_depth_nonzero, mean_depth_within2sd, median, proportion_of_reference_covered)
            total_depth <- depth
            total_presence <- presence
        } else {
            output <- rbind(output, cbind(bamfile, mean_depth, sd_depth, mean_depth_nonzero, mean_depth_within2sd, median, proportion_of_reference_covered))
            total_depth <- total_depth + depth
            total_presence <- total_presence + presence
        }
        }
        print(output)
        write_csv(output, path="/home/projects/dp_00007/people/hmon/EUostrea/02_data/Depth/output.WADD.csv")  #change path
        output2 <- output %>%
        mutate(across(where(is.numeric), round, 3))%>% 
        write_csv(output2, file = "/home/projects/dp_00007/people/hmon/EUostrea/02_data/Depth/samplespe_per_base_depth_presenceData.WADD.csv")
```
</details>



