#Load all required modules for the job
module load tools
module load ngs
module load jre/1.8.0-openjdk
module load trimmomatic/0.38

#Global variables
base=MOLU_03_EKDL200012859-1a-AK11200-AK15590_HNLVYDSXY_L1
ADAPTERS=home/projects/dp_00007/people/hmon/Shucking/01_infofiles/NexteraPE_NT.fa

cd /home/projects/dp_00007/people/hmon/Shucking
java -jar /services/tools/trimmomatic/0.38/trimmomatic-0.38.jar PE \
        -threads 8 \
        -phred33 \
	02_data/"$base"_1.fq.gz \
	02_data/"$base"_2.fq.gz \
	03_trimmed/"$base"_1.paired.fq.gz \
	03_trimmed/"$base"_1.single.fq.gz \
	03_trimmed/"$base"_2.paired.fq.gz \
	03_trimmed/"$base"_2.single.fq.gz \
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
base=MOLU_03_EKDL200012859-1a-AK11200-AK15590_HNLVYDSXY_L1

#Align reads
    echo "Aligning $base"
    ID=$(echo "@RG\tID:$base\tSM:$base\tPL:Illumina")

  #Align reads 1 step
    bwa mem -t "$NCPU" \
        -R "$ID" \
        "$GENOME" \
        03_trimmed/"$base"_1.paired.fq.gz 03_trimmed/"$base"_2.paired.fq.gz > 04_mapped/"$base".sam

        # Create bam file
    echo "Creating bam for $base"

    samtools view -bS -h -q 20 -F 4 \
    04_mapped/"$base".sam >04_mapped/"$base".bam


     echo "Creating sorted bam for $base"
        samtools sort 04_mapped/"$base".bam -o 04_mapped/"$base".sort.minq20.bam
        samtools index 04_mapped/"$base".sort.minq20.bam

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
I=04_mapped/"$base".sort.minq20.bam \
O=05_dedup/"$base".nocig.dedup.minq20.bam \
M=05_dedup/"$base".duprmmetrics.txt \
REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=SILENT

#scripts ClipOverlap with NO CIGAR on MarkDuplicates
/services/tools/bamutil/1.0.14/bam clipOverlap \
--in 05_dedup/"$base".nocig.dedup.minq20.bam \
--out 05_dedup/"$base".nocig.dedup_clipoverlap.minq20.bam \
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
samtools index 05_dedup/"$base".nocig.dedup_clipoverlap.minq20.bam

#Create list of potential in-dels nocig
java -jar /services/tools/gatk/3.8-0/GenomeAnalysisTK.jar \
-T RealignerTargetCreator \
-R $GENOME \
-I 05_dedup/"$base".nocig.dedup_clipoverlap.minq20.bam  \
-o 06_realigned/"$base".all_samples_for_indel_realigner.nocig.minq20.intervals

#Run the indel realigner tool nocig
java -jar /services/tools/gatk/3.8-0/GenomeAnalysisTK.jar \
-T IndelRealigner \
-R $GENOME \
-I 05_dedup/"$base".nocig.dedup_clipoverlap.minq20.bam \
-targetIntervals 06_realigned/"$base".all_samples_for_indel_realigner.nocig.minq20.intervals \
--consensusDeterminationModel USE_READS  --nWayOut _minq20.nocig.realigned.bam
