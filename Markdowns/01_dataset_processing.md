DATASET PROCESSING
================

- [Downsize MORL and USAM fastq files.](#Downsize-MORL-and-USAM-fastq-files.)
- [Processing USAM and MORL fastq.](#processing-usam-morl-fastq)
- [Fastq statistics.](#Fastq-stats)


## Downsize MORL and USAM fastq files.
    module load tools ngs  
    module load jdk/19 openjdk/19 java/1.8.0-openjdk jre/1.8.0-openjdk 
    module load bbmap/38.90

# samplerate=0.2 Randomly output only this fraction of reads; 1 means sampling is disabled.
    DIRFQ=/home/projects/dp_00007/people/hmon/Novaseq_MLX_USA
# FOR USAM population
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
    done
# FOR MORL population
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
    done


## Processing USAM MORL fastq
<details>

<summary> click for preprocess code </summary>

```#!/bin/bash
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

## Fastq stats
